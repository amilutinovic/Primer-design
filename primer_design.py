#!/usr/bin/env python
# coding: utf-8

# In[1]:


import random
import math
import statistics as stats
from collections import Counter


# # Simple prototype

# In[2]:


COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s): 
    return s.translate(COMP)[::-1]

def within(v, lo, hi): 
    return lo <= v <= hi
    
def clean_seq(seq):
    return "".join(ch for ch in seq.upper() if ch in "ACGT")

def _seed_match(a, b, seed=4):
    if len(a) < seed or len(b) < seed:
        return False
    return a[-seed:].upper() == revcomp(b.upper())[:seed]


# In[3]:


class Individual:
    """
    Individual of genetic algorithm representing a primer pair on a template
    """
    def __init__(self, template, f_start, f_len, r_start, r_len):
        self.f_start = f_start 
        self.f_len = f_len
        self.r_start = r_start
        self.r_len = r_len
        self.template = template 
        self.kcount_template = self.kmer_uniqueness(self.template, k=10)
        self.fitness = self.calc_fitness(kcount=self.kcount_template)

    # Sequence helpers
    def forward(self):
        return self.template[self.f_start:self.f_start+self.f_len]

    def reverse(self):
        frag = self.template[self.r_start:self.r_start+self.r_len]
        return revcomp(frag)

    def amplicon_len(self):
        return (self.r_start + self.r_len) - self.f_start

    def is_valid(self, len_win=(18,25), amp_bounds=(120,300)):
        n = len(self.template)
        if not (0 <= self.f_start < n and 0 <= self.r_start < n): 
            return False
        if self.f_start + self.f_len > n or self.r_start + self.r_len > n:
            return False
        if not (self.f_start < self.r_start): 
            return False
        if not within(self.f_len, *len_win): return False              
        if not within(self.r_len, *len_win): return False  
        if self.has_homopolymer(self.forward()) or self.has_homopolymer(self.reverse()): 
            return False
        amp = self.amplicon_len()
        if not within(amp, *amp_bounds): return False
        return True

    def repair(self, len_win=(18,25), amp_bounds=(120,300)):
        """
        Clamp lengths/indices into legal ranges and try to nudge r_start to satisfy amplicon bounds.
        """
        n = len(self.template)
        self.f_len = max(len_win[0], min(len_win[1], self.f_len))
        self.r_len = max(len_win[0], min(len_win[1], self.r_len))
        self.f_start = max(0, min(n - self.f_len, self.f_start))

        r_min = self.f_start + self.f_len + 1
        self.r_start = max(r_min, min(n - self.r_len, self.r_start))
        
        amp = self.amplicon_len()
        if amp < amp_bounds[0]:
            target = self.f_start + amp_bounds[0] - self.r_len
            self.r_start = max(r_min, min(n - self.r_len, target))
        elif amp > amp_bounds[1]:
            target = self.f_start + amp_bounds[1] - self.r_len
            self.r_start = max(r_min, min(n - self.r_len, target))
            
        
    # Scoring features
    def gc_content(self, s):
        g = s.count("G")+s.count("C")
        return 100.0 * g / max(1, len(s))

    def tm_wallace(self, s):
        a = s.count("A")+s.count("T")
        g = s.count("G")+s.count("C")
        return 2*a + 4*g

    # SantaLucia-lite nearest-neighbor Tm (more continuous than Wallace)
    @staticmethod 
    def _nn_dH_dS():
        # kcal/mol and cal/(mol*K) for DNA/DNA (SantaLucia 1998, simplified)
        # values are for dinucleotides in 5'->3' (nearest neighbor)
        return {
            "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
            "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
            "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
            "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
            "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
            "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
            "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
            "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
        }

    def tm_nn(self, s, Na=50e-3, primer_conc=500e-9):
        """
        Nearest-neighbor Tm (°C), monovalent salt correction (Owczarzy/SantaLucia-lite).
        Na: mol/L; primer_conc: total strand conc (mol/L) for non-self-complementary (Ct/4 term).
        """
        s = s.upper()
        if len(s) < 2:
            return 0.0
        dH = 0.0   # kcal/mol
        dS = 0.0   # cal/(mol*K)
        nn = self._nn_dH_dS()
        # initiation terms (rough)
        init_dH = 0.2
        init_dS = -5.7
        dH += init_dH
        dS += init_dS
        for i in range(len(s)-1):
            dimer = s[i:i+2]
            if dimer in nn:
                h, s_ = nn[dimer]
                dH += h
                dS += s_
        R = 1.987  # cal/(mol*K)
        # concentration term (Ct/4 for non-self-complementary)
        Ct = max(primer_conc, 1e-12)
        # base Tm in Kelvin
        TmK = (1000.0 * dH) / (dS + R * math.log(Ct/4.0))  # 1000 converts kcal->cal
        # monovalent salt correction
        TmK += (16.6 * math.log10(max(Na, 1e-6)))  # °C-ish in K units (approx)
        return TmK - 273.15


    def has_homopolymer(self, s, k=4):
        run = 1
        for i in range(1, len(s)):
            run = run+1 if s[i]==s[i-1] else 1
            if run >= k: return True
        return False

    def gc_clamp_score(self, s):
        tail = s[-5:] if len(s)>=5 else s
        gcs = sum(1 for c in tail if c in "GC")
        # bonus for 1-2 GC in 3' tail, penalty if too many
        if 1 <= gcs <= 2: 
            return 1.0
        if gcs >= 4: 
            return -1.0
        return 0.0

    def short_3prime_complement_score(self, p, q, seed=4):
        """
        Penalize if 3' seeds strongly complement (self- or cross-dimers)
        """
        rp = revcomp(p)
        rq = revcomp(q)
        s1 = p[-seed:]
        s2 = q[-seed:]
        return 1.0 if s1 and s2 and s1 == rq[:seed] else 0.0

    @staticmethod    
    def overlap_3p_score(p, q, max_k=6, min_k=3):  
        if not p or not q:                          
            return 0.0                               
        rp_head = revcomp(q)[:max_k]               
        tail = p[-max_k:]                           
        best = 0                                    
        for k in range(max_k, min_k-1, -1):         
            if tail[-k:] == rp_head[:k]:            
                best = k                             
                break                                
        if best == 0:                                
            return 0.0                               
        gc_tail = sum(c in "GC" for c in tail[-best:])  
        return (best - (min_k - 1)) + 0.25 * gc_tail 

    @staticmethod 
    def _hairpin_3p_score(p, max_k=6, min_k=3):     
        if not p:                                  
            return 0.0                               
        rc = revcomp(p)                             
        tail = p[-max_k:]                           
        best = 0                                    
        search = rc[2:-2] if len(rc) > 4 else rc    
        for k in range(max_k, min_k-1, -1):         
            if tail[-k:] in search:                 
                best = k                             
                break                                
        if best == 0:                                
            return 0.0                               
        gc_tail = sum(c in "GC" for c in tail[-best:])  
        return (best - (min_k - 1)) + 0.25 * gc_tail    

    def dimer_penalty(self, p, q):
        return (self.overlap_3p_score(p, p) \
                + self.overlap_3p_score(q, q) \
                + self.overlap_3p_score(p, q) \
                + self.overlap_3p_score(q, p))

    def hairpin_penalty(self, p, seed=4):
        """
        Penalize if 3' seed is complementary to an internal segment
        """
        return self._hairpin_3p_score(p) 

    def kmer_uniqueness(self, seq, k=10):
        kmers = [seq[i:i+k] for i in range(0, max(0,len(seq)-k+1))]
        return Counter(kmers)

    def primer_uniqueness_score(self, prim, kcount, k=10):
        if len(prim) < k:
            return 0.5
        ks = [prim[i:i+k] for i in range(len(prim)-k+1)]
        repeats = sum(max(kcount.get(km, 0) - 1, 0) for km in ks)
        return 1.0 / (1.0 + repeats)

    # Fitness
    def calc_fitness(self, len_win=(18,25), tm_win=(58,62), 
                     gc_win=(40,60),amplicon_win=(120,300), 
                     kcount=None):

        if not self.is_valid(len_win, amplicon_win):
            self.repair(len_win, amplicon_win)
            if not self.is_valid(len_win, amplicon_win):
                self.fitness = -1e9
                return self.fitness
        
        fwd = self.forward()
        rev = self.reverse()

        if min(len(fwd), len(rev)) < 16:
            self.fitness = -1e9
            return self.fitness

        amp_len = self.amplicon_len()
        if amp_len < amplicon_win[0] or amp_len > amplicon_win[1]:
            center = 0.5 * (amplicon_win[0] + amplicon_win[1])
            amp_score = -abs((amp_len - center) / 50.0)
        else:
            amp_score = 1.0


        f_tm = self.tm_wallace(fwd)
        r_tm = self.tm_wallace(rev)
        f_gc = self.gc_content(fwd)
        r_gc = self.gc_content(rev)

        def win_score(v, lo, hi, width=5.0):
            if within(v, lo, hi): 
                return 1.0
            edge = lo if v < lo else hi 
            return max(0.0, 1.0 - abs(v-edge)/width)

        tm_score = 0.5*win_score(f_tm, *tm_win) + 0.5*win_score(r_tm, *tm_win)
        gc_score = 0.5*win_score(f_gc, *gc_win) + 0.5*win_score(r_gc, *gc_win)
        len_score = 0.5*win_score(len(fwd), *len_win) + 0.5*win_score(len(rev), *len_win)
        clamp_score = 0.5*self.gc_clamp_score(fwd) + 0.5*self.gc_clamp_score(rev)

        # penalties
        run_pen = 1.0 if (self.has_homopolymer(fwd) or self.has_homopolymer(rev)) else 0.0
        dimer_pen = self.dimer_penalty(fwd, rev)
        hair_pen = self.hairpin_penalty(fwd) + self.hairpin_penalty(rev)

        if kcount is None:
            kcount = self.kcount_template
        uniq_score = 0.0
        if kcount is not None:
            uniq_score = 0.5* self.primer_uniqueness_score(fwd, kcount) + 0.5* self.primer_uniqueness_score(rev, kcount)

         # weights
        w_tm, w_gc, w_len, w_clamp, w_amp, w_uniq = 1.0, 1.0, 0.6, 0.5, 2.0, 0.8
        p_dimer, p_hair, p_run = 3.0, 2.0, 1.5

        fitness = (w_tm*tm_score + w_gc*gc_score + w_len*len_score \
                   +  w_clamp*clamp_score + w_amp*amp_score + w_uniq*uniq_score \
                   - p_dimer*dimer_pen - p_hair*hair_pen - p_run*run_pen)
        
        self.fitness = float(fitness)
        return self.fitness

    def genes(self):
        return [self.f_start, self.f_len, self.r_start, self.r_len]
        
    @staticmethod
    def from_genes(template, genes):
        return Individual(template, genes[0], genes[1], genes[2], genes[3])      

    def basic_metrics(self):
        F, R = self.forward(), self.reverse()
        return {
            "TmF": self.tm_nn(F), "TmR": self.tm_nn(R),
            "GCF": self.gc_content(F), "GCR": self.gc_content(R),
            "DG":  (self.hairpin_penalty(F) + self.hairpin_penalty(R) + self.dimer_penalty(F,R)),
            "Amp": self.amplicon_len(),
            "ClampF": self.gc_clamp_score(F), "ClampR": self.gc_clamp_score(R),
        }

        


# In[4]:


def selection(population, turnament_size):
    k = min(len(population), turnament_size)
    participants = random.sample(population, k)
    return max(participants, key=lambda x: x.fitness)


# In[5]:


def crossover(p1, p2, template, len_win=(18,25), amp_bounds=(120,300), try_cuts=6):
    g1 = p1.genes()
    g2 = p2.genes()
    
    M = min(len(g1), len(g2))
    cut = 1 if M == 1 else random.randint(1, M-1)
    
    for _ in range(max(1, try_cuts)):
        cut = random.randint(1, M-1)
        child_genes = g1[:cut] + g2[cut:M]
        child = Individual.from_genes(template, child_genes)
        child.repair(len_win, amp_bounds)
        if child.is_valid(len_win, amp_bounds):
            return child
            
    return p1 if p1.fitness >= p2.fitness else p2


# In[6]:


def mutation(child, len_win, amp_bounds, max_shift=5, rate=0.35):
    if random.random() < rate: 
        child.f_len += random.choice([-1,1])
    if random.random() < rate: 
        child.r_len += random.choice([-1,1])
    if random.random() < rate: 
        child.f_start += random.randint(-max_shift, max_shift)
    if random.random() < rate: 
        child.r_start += random.randint(-max_shift, max_shift)
    child.repair(len_win, amp_bounds)
    return child


# In[7]:


def random_individual(template, len_win, amp_bounds, attempts=400):
    template = clean_seq(template)
    n = len(template)
    for _ in range(attempts):
        f_len = random.randint(*len_win)
        r_len = random.randint(*len_win)
        f_start = random.randint(0, n - f_len - 1)
        min_gap = amp_bounds[0] - f_len
        r_min = max(f_start + f_len + 1, f_start + f_len + min_gap)
        r_min = min(r_min, n - r_len)
        if r_min > n - r_len: 
            continue
        r_start = random.randint(r_min, n - r_len)
        ind = Individual(template, f_start, f_len, r_start, r_len)
        if ind.is_valid(len_win, amp_bounds):
            ind.calc_fitness()
            return ind
    # Fallback
    ind = Individual(template, 10, 20, min(n-20, 10+20+150), 20)
    ind.repair(len_win, amp_bounds)
    ind.calc_fitness()
    return ind
    


# In[8]:


def ga(template,
       elitism=2,
       generations=120,
       pop_size=150,
       turnament_size=13,
       len_win=(18, 25),
       amp_bounds=(120, 300),
       seed=42):
    
    random.seed(seed)
    # Init population
    population = [random_individual(template, len_win, amp_bounds)for _ in range(pop_size)]

    for _ in range(generations):
        # Elitism
        population.sort(key=lambda ind: ind.fitness, reverse=True)
        next_population = population[:elitism]
        
        while len(next_population) < pop_size:
            p1 = selection(population, turnament_size=2)
            p2 = selection(population, turnament_size=2)
            child = crossover(p1, p2, template, len_win, amp_bounds)
            child = mutation(child, len_win=len_win, amp_bounds=amp_bounds, rate=0.35)

            next_population.append(child)
        population = next_population
    
    population.sort(key=lambda ind: ind.fitness, reverse=True)
    best = population[0]
    return best, population[:10]
    


# In[9]:


# template = ("ATGCGTACGTTGACCTGATCGATCGGATCCGATGCTAGCTAGCTAGGCTTACGATCGATCG"*6)
# best, top10 = ga(template, pop_size=120, generations=120, seed=42)
# print("BEST fitness:", round(best.fitness,3))
# print("BEST primers:")
# print("  F:", best.forward(), "| start:", best.f_start, "len:", best.f_len)
# print("  R:", best.reverse(), "| start:", best.r_start, "len:", best.r_len)
# print("  Amplicon:", best.amplicon_len(), "bp")


# # Multiplex NSGA-II for PCR primer design

# In[10]:


class MultiplexChromosome:
    def __init__(self, pairs, template):
        self.pairs = pairs
        self.template = template

    def clone(self):
        new_pairs = [Individual( self.template, p.f_start, p.f_len, p.r_start, p.r_len) for p in self.pairs]
        return MultiplexChromosome(new_pairs, self.template)

    def is_valid(self, len_win=(18,25), amp_bounds=(120,300), seed=4):
        for p in self.pairs:
            if not p.is_valid(len_win, amp_bounds):
                return False
        # global 3' seed conflicts (inter-pair)
        all_primers = [pp.forward() for pp in self.pairs] + [pp.reverse() for pp in self.pairs]
        for i in range(len(all_primers)):
            for j in range(i+1, len(all_primers)):
                if _seed_match(all_primers[i], all_primers[j], seed): 
                    return False
                if _seed_match(all_primers[j], all_primers[i], seed): 
                    return False
        return True

    def repair(self, len_win=(18,25), amp_bounds=(120,300)):
        for p in self.pairs: p.repair(len_win, amp_bounds)

    def evaluate(self, t_target=60.0, gc_target=50.0, amp_target=200, d_min=40.0):
        """Returns objectives tuple: (f_Tm, f_GC, f_DG, f_Amp)
        f_Tm: mean |Tm_pair - t_target| + (maxTm - minTm)
        f_GC: mean |GC - 50| across primers
        f_DG: sum of hairpin/self-dimer per pair + cross-dimers across ALL primers
        f_Amp: mean |Amp - amp_target| + penalty if closest amplicon gap < d_min
        """
        per = [p.basic_metrics() for p in self.pairs]
        pair_Tms = [ (m["TmF"]+m["TmR"])/2.0 for m in per ]
        mean_Tm_dev = sum(abs(t - t_target) for t in pair_Tms)/len(pair_Tms)
        span_Tm = stats.pstdev(pair_Tms) if len(pair_Tms) > 1 else 0.0
        f_Tm = mean_Tm_dev + span_Tm

        # GC deviation across all primers (both F and R)
        gc_vals = []
        for m in per: gc_vals.extend([m["GCF"], m["GCR"]])
        f_GC = sum(abs(g - gc_target) for g in gc_vals)/len(gc_vals)

        # ΔG-like penalties: intra-pair + global inter-primer cross
        intra = sum(m["DG"] for m in per)
        # all primers list
        all_prims = []
        for pp in self.pairs:
            all_prims.append(pp.forward()); all_prims.append(pp.reverse())
        inter = 0.0
        for i in range(len(all_prims)):
            for j in range(i+1, len(all_prims)):
                # cross-dimer check both directions
                inter += Individual.overlap_3p_score(all_prims[i], all_prims[j])
                inter += Individual.overlap_3p_score(all_prims[j], all_prims[i])
        f_DG = intra + inter

        # Amplicon objectives: closeness to target + spacing
        amps = [m["Amp"] for m in per]
        mean_amp_dev = sum(abs(a - amp_target) for a in amps)/len(amps)
        closest_gap = math.inf
        for i in range(len(amps)):
            for j in range(i+1, len(amps)):
                gap = abs(amps[i]-amps[j])
                if gap < closest_gap: closest_gap = gap
        penal = 0.0
        if len(amps) > 1 and closest_gap < d_min:
            penal = (d_min - closest_gap)  # linear penalty if too close
        f_Amp = mean_amp_dev + penal

        return (round(f_Tm,6), round(f_GC,6), round(f_DG,6), round(f_Amp,6))

    def genes(self):
        g=[]
        for p in self.pairs: g.extend(p.genes())
        return g 

    @staticmethod
    def from_genes_flat(flat, M, template):
        pairs=[]
        for i in range(M):
            block = flat[4*i:4*(i+1)]
            pairs.append(Individual.from_genes(template, block))
        return MultiplexChromosome(pairs, template)


# In[11]:


def multiplex_mutation(chromosome, rate=0.3, max_shift=5, len_win=(18,25), amp_bounds=(120,300)):
    # Pick one pair at random and mutate
    idx = random.randrange(len(chromosome.pairs))
    p = chromosome.pairs[idx]
    mutation(p, len_win, amp_bounds, max_shift, rate)


# In[12]:


def multiplex_crossover(p1, p2, len_win=(18,25), amp_bounds=(120,300), try_cuts=5):
    # Safety check
    if p1.template != p2.template:
        raise ValueError("Templates differ; cannot crossover safely.")
    if not p1.pairs or not p2.pairs:
        return p1.clone()

    M = min(len(p1.pairs), len(p2.pairs))
    
    # 1-point over pairs (swap sublists)
    for _ in range(max(1, try_cuts)):
        cut = 1 if M == 1 else random.randint(1, M-1)

        # IMPORTANT: pass (template, genes) to from_genes
        left  = [Individual.from_genes(p1.template, pp.genes()) for pp in p1.pairs[:cut]]
        right = [Individual.from_genes(p2.template, pp.genes()) for pp in p2.pairs[cut:M]]

        child = MultiplexChromosome(left + right, p1.template)
        child.repair(len_win=len_win, amp_bounds=amp_bounds)

        if child.is_valid(len_win=len_win, amp_bounds=amp_bounds):
            return child
            
    return p1.clone() if p1.evaluate() >= p2.evaluate() else p2.clone()


# In[13]:


def random_multiplex(template, M=3, len_win=(18,25), amp_bounds=(120,300)):
    pairs=[]
    for _ in range(M):
        pairs.append(random_individual(template, len_win, amp_bounds))
    mux = MultiplexChromosome(pairs, template)
    mux.repair(len_win, amp_bounds)
    return mux


# ### NSGA-II core

# In[14]:


def dominates(a,b):
    return all(x<=y for x,y in zip(a,b)) and any(x<y for x,y in zip(a,b))


# In[15]:


def fast_non_dominated_sort(F):
    S = [set() for _ in F] # S[p] all q that p dominates
    n = [0]*len(F) # Number of q that dominates over p
    fronts = [[]]
    for p in range(len(F)):
        for q in range(len(F)):
            if p==q: continue
            if dominates(F[p], F[q]): S[p].add(q)
            elif dominates(F[q], F[p]): n[p]+=1
        if n[p]==0: fronts[0].append(p)
    i=0
    while fronts[i]:
        Q=[]
        for p in fronts[i]:
            for q in S[p]:
                n[q]-=1
                if n[q]==0: Q.append(q)
        i+=1
        fronts.append(Q)
    if not fronts[-1]: fronts.pop()
    return fronts


# In[16]:


# If inside same Pareto front, choose a more diverse solution, keeping "the witdth" of Pareto front
def crowding_distance(F, idxs):
    if not idxs: return {}
    m = len(F[0])
    dist = {i:0.0 for i in idxs}
    for k in range(m):
        s = sorted(idxs, key=lambda i: F[i][k])
        fmin, fmax = F[s[0]][k], F[s[-1]][k]
        dist[s[0]] = dist[s[-1]] = float("inf") # Distance from neighbours
        if fmax==fmin: continue
        for r in range(1,len(s)-1):
            i_prev, i_next, i = s[r-1], s[r+1], s[r]
            dist[i] += (F[i_next][k]-F[i_prev][k])/(fmax-fmin)
    return dist


# In[17]:


# Selecting parents for next generation, lower rank is better Pareto front, if == then crowding distance
def crowding_tournament(a, b, rank, crowding_distance):
    if rank[a] < rank[b]: return a
    if rank[b] < rank[a]: return b
    return a if crowding_distance.get(a,0.0) >= crowding_distance.get(b,0.0) else b


# In[18]:


def nsga2_multiplex(template, M=3, pop_size=120, generations=120,
                    len_win=(18,25), amp_bounds=(120,300),
                    t_target=60.0, gc_target=50.0, amp_target=None, d_min=40.0, seed=42):
    random.seed(seed)
    template = clean_seq(template)   
    if amp_target is None:
        amp_target = (amp_bounds[0]+amp_bounds[1])//2

    # Init population
    P = [random_multiplex(template, M, len_win, amp_bounds) for _ in range(pop_size)]

    # Evaluate initial population
    def eval_population(population):
        objs=[]
        for multi_ch in population:
            if not multi_ch.is_valid(len_win, amp_bounds):
                objs.append((math.inf, math.inf, math.inf, math.inf))
            else:
                objs.append(multi_ch.evaluate(t_target, gc_target, amp_target, d_min))
        return objs

    F = eval_population(P)

    for _ in range(generations):
        fronts = fast_non_dominated_sort(F)
        rank={}
        for r, fr in enumerate(fronts):
            for i in fr: rank[i]=r
        crowd={}
        for fr in fronts:
            crowd.update(crowding_distance(F, fr))

        # Selection
        parents=[]
        while len(parents) < pop_size:
            i,j = random.randrange(pop_size), random.randrange(pop_size)
            parents.append(P[crowding_tournament(i,j,rank,crowd)])

        # variation
        children =[]
        for i in range(0, pop_size, 2):
            p1 = parents[i]
            p2 = parents[(i+1)%pop_size]
            child1 = multiplex_crossover(p1, p2, len_win, amp_bounds)
            multiplex_mutation(child1, rate=0.35, len_win=len_win, amp_bounds=amp_bounds)
            child2 = multiplex_crossover(p2, p1, len_win, amp_bounds)
            multiplex_mutation(child2, rate=0.35, len_win=len_win, amp_bounds=amp_bounds)
            children.extend([child1, child2])
        children = children[:pop_size]

        # Elitist replacement
        R = P + children # popula + children 
        FR = eval_population(R)
        frontsR = fast_non_dominated_sort(FR)
        P_next=[]; F_next=[]
        for fr in frontsR:
            if len(P_next)+len(fr) <= pop_size:
                # Adding the whole front if there is room
                P_next.extend([R[i] for i in fr])
                F_next.extend([FR[i] for i in fr])
            else:
                # There is no room -> crowding “truncation”
                cd = crowding_distance(FR, fr)
                fr_sorted = sorted(fr, key=lambda i: cd[i], reverse=True)
                slots = pop_size - len(P_next)
                take = fr_sorted[:slots]
                P_next.extend([R[i] for i in take])
                F_next.extend([FR[i] for i in take])
                break
        P, F = P_next, F_next

    # Final Pareto
    fronts = fast_non_dominated_sort(F)
    pareto_idx = fronts[0] if fronts else []
    pareto = []
    for i in pareto_idx:
        ch = P[i]
        objs = F[i]  # Tuple: (f_Tm, f_GC, f_DG, f_Amp)

        row = {"f_Tm": objs[0], "f_GC": objs[1], "f_DG": objs[2], "f_Amp": objs[3]}
        
        for pi, pp in enumerate(ch.pairs, start=1):
            fwd = pp.forward()
            rev = pp.reverse()
            row.update({
                f"p{pi}_f_start": pp.f_start,  f"p{pi}_f_len": pp.f_len,
                f"p{pi}_r_start": pp.r_start,  f"p{pi}_r_len": pp.r_len,
                f"p{pi}_F": fwd,               f"p{pi}_R": rev,
                f"p{pi}_Amp": pp.amplicon_len(),  
                f"p{pi}_TmF": pp.tm_wallace(fwd), f"p{pi}_TmR": pp.tm_wallace(rev),
                f"p{pi}_GCF": round(pp.gc_content(fwd), 1),
                f"p{pi}_GCR": round(pp.gc_content(rev), 1),
                f"p{pi}_ClampF": pp.gc_clamp_score(fwd),
                f"p{pi}_ClampR": pp.gc_clamp_score(rev),
            })
        pareto.append(row)
    # Sort by Tm harmony then DG
    pareto.sort(key=lambda r: (r["f_Tm"], r["f_DG"], r["f_Amp"], r["f_GC"]))
    return pareto


# In[19]:


# Treats two rows as the same if all four objectives (rounded) and each pair’s (F, R, Amp) match.
def dedup_pareto(rows, M=3, nd=6):
    seen = set()
    uniq = []
    for r in rows:
        sig = (
            round(r["f_Tm"], nd), round(r["f_GC"], nd),
            round(r["f_DG"], nd), round(r["f_Amp"], nd)
        ) + tuple((r[f"p{i}_F"], r[f"p{i}_R"], r[f"p{i}_Amp"]) for i in range(1, M+1))
        if sig in seen:
            continue
        seen.add(sig)
        uniq.append(r)
    return uniq


# In[20]:


# template = ("ATGCGTACGTTGACCTGATCGATCGGATCCGATGCTAGCTAGCTAGGCTTACGATCGATCG"*8)

# pareto = nsga2_multiplex(
#        template=template,
#        M=3,
#        pop_size=150,
#        generations=150,
#        len_win=(18,25),
#        amp_bounds=(120,300),
#        t_target=60.0,
#        gc_target=50.0,
#        amp_target=200,
#        d_min=40.0,
#        seed=42
#    )

# pareto = dedup_pareto(pareto, M=3)

# for i, row in enumerate(pareto[:3], start=1):
    # print(f"\n=== Pareto #{i} ===")
    # print({k: row[k] for k in ["f_Tm","f_GC","f_DG","f_Amp"]})
    # for p in range(1,4):
        # print(f" Pair {p}: Amp={row[f'p{p}_Amp']}  F={row[f'p{p}_F']}  R={row[f'p{p}_R']}")


# ## Pareto #1
# - f_Tm = 3.883 → composite Tm cost ≈ 3.9: pairs are a few °C off 60 °C and/or not perfectly harmonized across pairs.
# - f_GC = 5.694 → average deviation from 50% GC is ~5.69 percentage points (moderate; can be improved).
# - f_DG = 0.0 → no 3′-overlap/self-dimer/hairpin signal by the current graded check (good safety).
# - f_Amp = 26.67 → amplicons 156/147/206 with a single target of 200 bp: mean |Δ| = (|156−200|+|147−200|+|206−200|)/3 = (44+53+6)/3 = 34.33; closest gap = min( |156−147|=9, |156−206|=50, |147−206|=59 ) = 9 bp < d_min(40) ⇒ spacing penalty 31;
# total = 34.33 + 31 = 65.33.
#   
# Verdict: Dimer-safe but spacing is poor (147 vs 156 only 9 bp apart) and sizes are far from 200 for two pairs.
# 
# ## Pareto #2
# 
# - f_Tm = 3.959 → similar overall Tm mistuning/spread as #1.
# - f_GC = 5.694 → average GC off by ~5.69 pp (same as #1).
# - f_DG = 7.5 → elevated dimer/hairpin risk; multiple 3′ overlaps contribute (needs sequence/end adjustments).
# - f_Amp = 46.0 → amplicons 152/173/206 with target 200 bp: mean |Δ| = (|152−200|+|173−200|+|206−200|)/3 = (48+27+6)/3 = 27.0;
# closest gap = min( |152−173|=21, |152−206|=54, |173−206|=33 ) = 21 bp < 40 ⇒ spacing penalty 19; total = 27 + 19 = 46.0.
#   
# Verdict: Best of the three on amplicon targeting/spacing, but unsafe (high f_DG).
# 
# ## Pareto #3
# 
# - f_Tm = 4.081 → slightly worse Tm composite than #1/#2 (a bit more off/spread).
# - f_GC = 5.028 → best GC of the three (≈5.03 pp from 50% on average).
# - f_DG = 0.0 → clean for dimers/hairpins by current check (good).
# - f_Amp = 60.0 → amplicons 160/147/206 with target 200 bp: mean |Δ| = (|160−200|+|147−200|+|206−200|)/3 = (40+53+6)/3 = 33.0; closest gap = min( |160−147|=13, |160−206|=46, |147−206|=59 ) = 13 bp < 40 ⇒ spacing penalty 27; total = 33 + 27 = 60.0.
# 
# Verdict: Dimer-safe and best GC balance, but spacing again too tight (147 vs 160 = 13 bp). Good starting point if you value safety;

# In[1]:


from Bio import SeqIO

def load_template(file_path):
    """Load template from FASTA or TXT file, or return raw string if not a path."""
    if file_path.endswith(".fasta"):
        record = next(SeqIO.parse(file_path, "fasta"))
        return str(record.seq)
    elif file_path.endswith(".txt"):
        with open(file_path) as f:
            return f.read().strip()
    else:
        return file_path


# In[3]:


def run_algorithm(
    template_file,
    mode="multiplex",
    M=3,
    pop_size=150,
    generations=150,
    len_win=(18,25),
    amp_bounds=(120,300),
    t_target=60.0,
    gc_target=50.0,
    amp_target=200,
    d_min=40.0,
    seed=42
):
    """Run GA algorithm and return results as text."""

    # load template
    template = load_template(template_file)


    pareto = nsga2_multiplex(
        template=template,
        M=M,
        pop_size=pop_size,
        generations=generations,
        len_win=len_win,
        amp_bounds=amp_bounds,
        t_target=t_target,
        gc_target=gc_target,
        amp_target=amp_target,
        d_min=d_min,
        seed=seed
    )
    pareto = dedup_pareto(pareto, M=M)

    output_lines = []
    for i, row in enumerate(pareto[:3], start=1):
        output_lines.append(f"\n=== Pareto #{i} ===")
        output_lines.append(str({k: row[k] for k in ["f_Tm","f_GC","f_DG","f_Amp"]}))
        for p in range(1, M+1):
            output_lines.append(
                f" Pair {p}: Amp={row[f'p{p}_Amp']}  F={row[f'p{p}_F']}  R={row[f'p{p}_R']}"
            )

    summary_str = "\n".join(output_lines)
    return summary_str


# In[ ]:





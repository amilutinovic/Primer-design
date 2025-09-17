import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QFileDialog, QMessageBox,
    QGridLayout, QGroupBox, QTextEdit
)
from primer_design import run_algorithm  # exported notebook as .py

class PrimerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PCR Primer Design")
        self.resize(700, 600)

        main_layout = QVBoxLayout(self)

        # template file
        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Template file:"))
        self.template_path = QLineEdit()
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_file)
        file_layout.addWidget(self.template_path)
        file_layout.addWidget(browse_btn)
        main_layout.addLayout(file_layout)

        # parameters
        params_group = QGroupBox("Parameters")
        grid = QGridLayout()
        params_group.setLayout(grid)

        defaults = {
            "M": 3, "pop_size": 150, "generations": 150,
            "len_win min": 18, "len_win max": 25,
            "amp_bounds min": 120, "amp_bounds max": 300,
            "t_target": 60.0, "gc_target": 50.0,
            "amp_target": 200, "d_min": 40.0
        }

        self.entries = {}
        row, col = 0, 0
        for key, val in defaults.items():
            lbl = QLabel(key)
            entry = QLineEdit(str(val))
            self.entries[key] = entry
            grid.addWidget(lbl, row, col)
            grid.addWidget(entry, row, col + 1)

            if col == 0:
                col = 2
            else:
                col = 0
                row += 1

        main_layout.addWidget(params_group)

        # run button
        self.run_btn = QPushButton("Run Algorithm")
        self.run_btn.clicked.connect(self.run_program)
        main_layout.addWidget(self.run_btn)

        # output area
        self.output_area = QTextEdit()
        self.output_area.setReadOnly(True)
        main_layout.addWidget(QLabel("Algorithm Output:"))
        main_layout.addWidget(self.output_area)

    def browse_file(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self, "Select Template File", "", "FASTA/TXT Files (*.fasta *.txt)"
        )
        if file_name:
            self.template_path.setText(file_name)

    def run_program(self):
        if not self.template_path.text():
            QMessageBox.warning(self, "Error", "Please select a template file!")
            return

        self.output_area.setPlainText("Running algorithm... Please wait.")
        QApplication.processEvents()  

        try:
            params = {
                "template_file": self.template_path.text(),
                "M": int(self.entries["M"].text()),
                "pop_size": int(self.entries["pop_size"].text()),
                "generations": int(self.entries["generations"].text()),
                "len_win": (
                    int(self.entries["len_win min"].text()),
                    int(self.entries["len_win max"].text())
                ),
                "amp_bounds": (
                    int(self.entries["amp_bounds min"].text()),
                    int(self.entries["amp_bounds max"].text())
                ),
                "t_target": float(self.entries["t_target"].text()),
                "gc_target": float(self.entries["gc_target"].text()),
                "amp_target": int(self.entries["amp_target"].text()),
                "d_min": float(self.entries["d_min"].text()),
            }

            
            result = run_algorithm(**params)
            self.output_area.setPlainText(result)

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PrimerGUI()
    window.show()
    sys.exit(app.exec_())

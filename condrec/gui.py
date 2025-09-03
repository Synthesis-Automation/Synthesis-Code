import json
import sys
import time
from typing import Optional

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QApplication,
        QMainWindow,
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QPlainTextEdit,
        QLineEdit,
        QPushButton,
        QFileDialog,
        QCheckBox,
        QLabel,
        QMessageBox,
        QSizePolicy,
    )
except Exception as e:  # pragma: no cover - only triggered when GUI deps missing
    qt_import_error = e
    QApplication = None  # type: ignore
else:
    qt_import_error = None

try:
    # When executed as a module: python -m condrec.gui
    from .analysis import analyze_reaction_smiles, format_analysis_summary
except Exception:
    # When executed directly: python condrec/gui.py
    try:
        import os
        pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        if pkg_root not in sys.path:
            sys.path.insert(0, pkg_root)
        from condrec.analysis import analyze_reaction_smiles, format_analysis_summary  # type: ignore
    except Exception as _e:  # pragma: no cover
        # Defer raising until main() so GUI error dialog/print can show guidance
        analyze_import_error = _e
    else:
        analyze_import_error = None
else:
    analyze_import_error = None


APP_TITLE = "CondRec – Reaction Analyzer"


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_TITLE)
        self.resize(900, 600)

        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # Display area (read-only)
        self.display = QPlainTextEdit(self)
        self.display.setReadOnly(True)
        self.display.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.display.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.display, stretch=1)

        # Bottom input row
        bottom = QHBoxLayout()
        layout.addLayout(bottom)

        self.json_checkbox = QCheckBox("JSON", self)
        self.json_checkbox.setToolTip("Show analysis as JSON instead of text")
        bottom.addWidget(self.json_checkbox)

        bottom.addWidget(QLabel("Reaction SMILES:", self))

        self.input = QLineEdit(self)
        self.input.setPlaceholderText("reactants>agents>products (agents optional)")
        self.input.returnPressed.connect(self.on_send)
        bottom.addWidget(self.input, stretch=1)

        self.open_btn = QPushButton("Open…", self)
        self.open_btn.clicked.connect(self.on_open_file)
        bottom.addWidget(self.open_btn)

        self.send_btn = QPushButton("Analyze", self)
        self.send_btn.setDefault(True)
        self.send_btn.clicked.connect(self.on_send)
        bottom.addWidget(self.send_btn)

        self.clear_btn = QPushButton("Clear", self)
        self.clear_btn.clicked.connect(self.display.clear)
        bottom.addWidget(self.clear_btn)

        self.statusBar().showMessage("Ready")

        self._append_welcome()

    def _append(self, text: str) -> None:
        self.display.appendPlainText(text)
        self.display.appendPlainText("")
        # Auto-scroll to bottom
        cursor = self.display.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        self.display.setTextCursor(cursor)

    def _append_welcome(self) -> None:
        self._append(
            (
                "Welcome! Enter a reaction SMILES below and press Analyze.\n"
                "Format: reactants>agents>products (agents optional).\n"
                "Tip: Use the JSON checkbox for machine-readable output."
            )
        )

    def on_open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open reaction SMILES file", "", "Text Files (*.txt);;All Files (*)")
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                content = f.read().strip()
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Failed to read file:\n{exc}")
            return
        self.input.setText(content)
        self.input.setFocus()

    def on_send(self) -> None:
        rxn = self.input.text().strip()
        if not rxn:
            self.statusBar().showMessage("Please enter a reaction SMILES.", 3000)
            return

        self._append(f">>> {rxn}")
        t0 = time.time()
        try:
            result = analyze_reaction_smiles(rxn)
        except ImportError as e:
            self._append(str(e))
            self.statusBar().showMessage("RDKit not available", 5000)
            return
        except Exception as e:
            self._append(f"Error: {e}")
            self.statusBar().showMessage("Failed to analyze", 5000)
            return

        elapsed_ms = int((time.time() - t0) * 1000)
        if self.json_checkbox.isChecked():
            out = json.dumps(result, indent=2, sort_keys=True)
        else:
            out = format_analysis_summary(result)
        self._append(out)
        self.statusBar().showMessage(f"Done in {elapsed_ms} ms", 3000)
        self.input.selectAll()


def main(argv: Optional[list[str]] = None) -> int:
    if qt_import_error is not None:
        print(
            "PyQt6 is required for the GUI. Install via: 'pip install PyQt6' or 'pip install .[gui]'.\n"
            f"Import error: {qt_import_error}",
            file=sys.stderr,
        )
        return 2

    if analyze_import_error is not None:
        print(
            "Failed to import analysis module. If you are running the file directly, try one of: \n"
            "  - 'python -m condrec.gui' (preferred)\n"
            "  - 'condrec-gui' (after 'pip install -e .')\n"
            f"Import error: {analyze_import_error}",
            file=sys.stderr,
        )
        return 2

    app = QApplication(argv or sys.argv)
    win = MainWindow()
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())

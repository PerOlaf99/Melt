"""Input dialogs for the varmelt GUI: paste sequences / wt-mut pairs."""
from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QFormLayout,
                               QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit,
                               QVBoxLayout)

from .. import webapp
from .model import Item, CLAMP_SIDES


def _ok_buttons(dialog, text="Add"):
    box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    box.button(QDialogButtonBox.Ok).setText(text)
    box.accepted.connect(dialog.accept)
    box.rejected.connect(dialog.reject)
    return box


class _ClampRow(QHBoxLayout):
    """A clamp selector ({5',3',none}) pre-filled with *default_*."""

    def __init__(self, default_side="5'"):
        super().__init__()
        from PySide6.QtWidgets import QComboBox
        self.combo = QComboBox()
        for s in CLAMP_SIDES:
            self.combo.addItem("GC clamp: 5' (left)" if s == "5'"
                               else "GC clamp: 3' (right)" if s == "3'"
                               else "GC clamp: none", s)
        self.combo.setCurrentIndex(CLAMP_SIDES.index(default_side))
        self.addWidget(self.combo)
        self.addStretch(1)

    def value(self):
        return self.combo.currentData()


class AddSequenceDialog(QDialog):
    """One or more amplicon strings (plain lines or FASTA); wt/mut name
    pairs become :class:`Item` pairs, everything else single sequences."""

    def __init__(self, parent=None, default_clamp="5'"):
        super().__init__(parent)
        self.setWindowTitle("Add DNA sequence(s)")
        self.setMinimumWidth(560)
        form = QFormLayout(self)
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("optional — headers are used per "
                                          "record when left empty")
        self.seq_edit = QPlainTextEdit()
        self.seq_edit.setPlaceholderText(
            ">frag1\nACGTTG…\n>frag2_wt\n…\n>frag2_mut\n…")
        self.seq_edit.setMinimumHeight(180)
        self.clamp_row = _ClampRow(default_clamp)
        help_lbl = QLabel("<span style='font-size:11px;color:#888'>"
                          "First 20 bases = forward primer; reverse "
                          "complement of last 20 = reverse primer.  Names "
                          "ending <b>_wt</b>/<b>_mut</b> (or _ref/_alt) with "
                          "the same base name are paired and plotted "
                          "together.</span>")
        help_lbl.setWordWrap(True)
        form.addRow("Name", self.name_edit)
        form.addRow("Sequence(s)", self.seq_edit)
        form.addRow("Clamp", self.clamp_row)
        form.addRow("", help_lbl)
        form.addRow("", _ok_buttons(self, "Add"))

    def items(self) -> list:
        """Parse the pasted text into GUI items (returned on accept)."""
        name = self.name_edit.text().strip()
        text = self.seq_edit.toPlainText()
        records = webapp._parse_fasta(text)
        out = []
        for item in webapp._pair_sequences(records):
            if item["kind"] == "seq":
                out.append(Item(kind="seq", name=item["name"] or name or
                                "sequence", seq=item["seq"],
                                clamp=self.clamp_row.value()))
            else:
                out.append(Item(kind="pair", name=item["name"] or name or
                                "pair", wt=item["wt"], mut=item["mut"],
                                clamp=self.clamp_row.value()))
        return out


class AddPairDialog(QDialog):
    """Explicit wildtype + mutant pair."""

    def __init__(self, parent=None, default_clamp="5'"):
        super().__init__(parent)
        self.setWindowTitle("Add wildtype / mutant pair")
        self.setMinimumWidth(560)
        form = QFormLayout(self)
        self.name_edit = QLineEdit()
        self.wt_edit = QPlainTextEdit()
        self.wt_edit.setMinimumHeight(110)
        self.mut_edit = QPlainTextEdit()
        self.mut_edit.setMinimumHeight(110)
        self.clamp_row = _ClampRow(default_clamp)
        form.addRow("Name", self.name_edit)
        form.addRow("Wildtype (5'->3')", self.wt_edit)
        form.addRow("Mutant (5'->3')", self.mut_edit)
        form.addRow("Clamp", self.clamp_row)
        form.addRow("", _ok_buttons(self, "Add"))

    def item(self) -> Item:
        wt = self.wt_edit.toPlainText()
        mut = self.mut_edit.toPlainText()
        return Item(kind="pair", name=self.name_edit.text().strip() or "pair",
                    wt=wt, mut=mut, clamp=self.clamp_row.value())
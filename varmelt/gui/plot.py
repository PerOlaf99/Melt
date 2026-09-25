"""QPainter-based melting-profile chart for the varmelt GUI.

Draws the per-base local melting temperature (WinMelt-style) as polylines:
black = wildtype/reference, dashed red = mutant/variant.  Also shades the
GC-clamp oligo on the end it occupies and marks the mutation base.
"""
import math

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QSizePolicy, QWidget


class MeltChart(QWidget):
    """Chart widget; call :meth:`set_map` to (re)plot an item's map."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(620, 260)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setAutoFillBackground(True)
        pal = self.palette()
        pal.setColor(self.backgroundRole(), QColor("white"))
        self.setPalette(pal)
        self._seq = ""
        self._ref = []
        self._alt = []
        self._clamp_len = 0
        self._clamp_side = None
        self._mark_idx = None
        self._mark_label = "mutation"
        self._indel = 0
        self._xmin, self._xmax = 30.0, 98.0
        self._paired = False
        self._status = "no data"

    def clear(self, status="no data"):
        self._seq = ""
        self._ref = []
        self._alt = []
        self._clamp_len = 0
        self._clamp_side = None
        self._mark_idx = None
        self._paired = False
        self._status = status
        self.update()

    def set_map(self, ref_prof, alt_prof=None, seq="", clamp_len=0,
                clamp_side=None, mark_idx=None, mark_label="mutation",
                indel=0, paired=False):
        self._ref = list(ref_prof) if ref_prof is not None else []
        self._alt = (list(alt_prof)
                     if alt_prof is not None
                     and list(alt_prof) != list(ref_prof) else [])
        self._seq = seq
        self._clamp_len = max(0, int(clamp_len))
        self._clamp_side = clamp_side
        self._mark_idx = mark_idx
        self._mark_label = mark_label
        self._indel = int(indel)
        self._paired = bool(paired)
        self._status = "ready"
        self.update()

    # -- geometry / helpers ----------------------------------------------- #
    def _xf(self, g, pad_l, plot_w, n):
        return pad_l + g * plot_w / max(n - 1, 1)

    def _yf(self, t, pad_t, plot_h):
        return pad_t + (1.0 - (max(float(t), self._xmin) - self._xmin)
                        / (self._xmax - self._xmin)) * plot_h

    def _draw_path(self, qp, prof, g, n, pad_l, pad_t, plot_w, plot_h, width):
        pen = qp.pen()
        pen.setWidthF(width)
        qp.setPen(pen)
        seg = []
        for i, t in prof:
            if math.isnan(t):
                if len(seg) >= 2:
                    qp.drawPolyline(QPolygonF(seg))
                seg = []
                continue
            seg.append(QPointF(self._xf(g(i), pad_l, plot_w, n),
                               self._yf(t, pad_t, plot_h)))
        if len(seg) >= 2:
            qp.drawPolyline(QPolygonF(seg))

    # -- painting --------------------------------------------------------- #
    def paintEvent(self, event):                        # noqa: N802
        qp = QPainter(self)
        qp.setRenderHint(QPainter.Antialiasing)
        w, h = self.width(), self.height()
        qp.fillRect(0, 0, w, h, QColor("white"))
        pad_l, pad_r, pad_t, pad_b = 46, 12, 10, 26
        plot_w, plot_h = w - pad_l - pad_r, h - pad_t - pad_b

        if not self._ref:
            qp.setPen(QColor("#999"))
            qp.drawText(self.rect(), Qt.AlignCenter, self._status or "no data")
            return

        n = len(self._ref)

        # grid + y labels
        qp.setFont(QFont("Helvetica", 8))
        for t in range(40, 96, 5):
            y = self._yf(t, pad_t, plot_h)
            qp.setPen(QColor("#ddd"))
            qp.drawLine(int(pad_l), int(y), int(w - pad_r), int(y))
            qp.setPen(QColor("#777"))
            qp.drawText(int(pad_l - 42), int(y - 7), 38, 14,
                        Qt.AlignRight | Qt.AlignVCenter, f"{t} °C")

        # clamp oligo shade on the end it actually occupies
        if self._clamp_len and self._clamp_len > 0:
            cl = min(self._clamp_len, n)
            x0 = pad_l
            x1 = pad_l + plot_w
            if self._clamp_side == "3'":
                x0 = pad_l + (n - cl) * plot_w / max(n - 1, 1)
                x1 = pad_l + plot_w
            else:
                x1 = pad_l + cl * plot_w / max(n - 1, 1)
            qp.setPen(Qt.NoPen)
            qp.setBrush(QColor("#f0eded"))
            qp.drawRect(int(x0), int(pad_t), int(max(x1 - x0, 2)),
                        int(plot_h))
            qp.setPen(QColor("#a99"))
            qp.drawText(int((x0 + x1) / 2 - 45), int(pad_t + 4), 90, 14,
                        Qt.AlignCenter, "GC clamp")

        def g_alt(i):
            v = self._mark_idx if self._mark_idx is not None else 0
            return i if i < v else i + self._indel

        # wildtype (black)
        qp.setBrush(Qt.NoBrush)
        qp.setPen(QColor("#222"))
        self._draw_path(qp, ((i, t) for i, t in enumerate(self._ref)),
                        lambda i: i, n, pad_l, pad_t, plot_w, plot_h, 1.6)

        # mutant (dashed red) unless it coincides with the wildtype
        if self._alt:
            pen = QPen(QColor("#C22"))
            pen.setWidthF(1.2)
            pen.setStyle(Qt.DashLine)
            qp.setPen(pen)
            self._draw_path(qp, ((i, t) for i, t in enumerate(self._alt)),
                            g_alt, n, pad_l, pad_t, plot_w, plot_h, 1.2)

        # mutation / variant mark
        if self._mark_idx is not None and 0 <= self._mark_idx < n:
            x = self._xf(self._mark_idx, pad_l, plot_w, n)
            qp.setPen(QPen(QColor("#D22")))
            qp.drawLine(int(x), int(pad_t), int(x), int(pad_t + plot_h))
            qp.setPen(QColor("#D22"))
            qp.setFont(QFont("Helvetica", 8, QFont.Bold))
            qp.drawText(int(x + 5), int(pad_t + 12), 90, 14,
                        Qt.AlignLeft | Qt.AlignVCenter,
                        self._mark_label or "variant")

        # footer caption
        qp.setFont(QFont("Helvetica", 8))
        qp.setPen(QColor("#999"))
        if self._paired and self._alt:
            cap = "black = wildtype    dashed red = mutant    (5'->3')"
        elif self._alt:
            cap = "black = reference    dashed red = variant    (5'->3')"
        else:
            cap = "reference melt map    (5'->3')"
        qp.drawText(int(pad_l), h - 20, plot_w, 16, Qt.AlignLeft, cap)
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
        self._delta_drawn = False

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
        self._delta_drawn = False
        self.update()

    # -- geometry / helpers ----------------------------------------------- #
    def _xf(self, g, pad_l, plot_w, n):
        return pad_l + g * plot_w / max(n - 1, 1)

    def _yf(self, t, pad_t, plot_h):
        return pad_t + (1.0 - (max(float(t), self._xmin) - self._xmin)
                        / (self._xmax - self._xmin)) * plot_h

    def _amp(self):
        """0-based fragment coordinates: the GC-clamp oligo is *not* part of
        the primer, so the primers live inside the amplicon which sits after
        the clamp on the 5' side and before it on the 3' side."""
        n = len(self._ref)
        cl = min(self._clamp_len, n) if self._clamp_side in ("5'", "3'") \
            else 0
        amp = n - cl
        offset = cl if self._clamp_side == "5'" else 0
        return offset, amp, n

    def primer_regions(self):
        """((fp0, fp1), (rp0, rp1)): half-open primer base intervals in
        fragment coordinates.  With a left GC clamp the forward primer
        starts at base ``len(GC_CLAMP)`` (0-based), i.e. base 43 1-based."""
        offset, amp, n = self._amp()
        return ((offset, min(offset + 20, n)),
                (max(offset + amp - 20, 0), offset + amp))

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

    def _draw_delta_band(self, qp, n, pad_l, pad_t, plot_w, plot_h):
        """Translucent area between the wt and mutant curves (aligned at the
        mutation, leaving a gap for an indel) -- the visual melt separation
        for a wt/mut pair.  Sets :attr:`_delta_drawn` when a band was
        actually rendered."""
        if not self._alt or not self._paired:
            return
        n_alt = len(self._alt)
        var = self._mark_idx if self._mark_idx is not None else 0

        def g_alt(j):
            return j if j < var else j + self._indel

        def flush(top, bot):
            if len(top) >= 2 and len(bot) >= 2:
                qp.setPen(Qt.NoPen)
                qp.setBrush(QColor(255, 195, 195, 150))
                qp.drawPolygon(QPolygonF(
                    [QPointF(x, y) for x, y in top]
                    + [QPointF(x, y) for x, y in reversed(bot)]))
                self._delta_drawn = True

        top, bot = [], []
        for i in range(n):
            rt = self._ref[i]
            j = i if i < var else i - self._indel
            if math.isnan(rt) or not (0 <= j < n_alt) \
                    or math.isnan(self._alt[j]):
                flush(top, bot)
                top, bot = [], []
                continue
            x = self._xf(i, pad_l, plot_w, n)
            yr = self._yf(rt, pad_t, plot_h)
            ya = self._yf(self._alt[j], pad_t, plot_h)
            if yr <= ya:
                top.append((x, yr))
                bot.append((x, ya))
            else:
                top.append((x, ya))
                bot.append((x, yr))
        flush(top, bot)

    def _draw_ruler(self, qp, n, pad_l, pad_r, plot_w):
        """Base numbers above the map; primer starts get a marked caret."""
        y0, y1 = 14, 20
        qp.setFont(QFont("Helvetica", 7))
        every = 20
        for b in range(0, n, every):
            x = self._xf(b, pad_l, plot_w, n)
            qp.setPen(QColor("#bbb"))
            qp.drawLine(int(x), y0, int(x), y1)
            qp.setPen(QColor("#888"))
            if n > 120 and b % (every * 2) and b != 0:
                continue                        # fewer labels when crowded
            qp.drawText(int(x - 16), y1 + 1, 32, 10, Qt.AlignHCenter,
                        str(b + 1))
        fp_reg, rp_reg = self.primer_regions()
        for x0, x1 in (fp_reg, rp_reg):
            if x1 <= x0:
                continue
            x = self._xf(x0, pad_l, plot_w, n)
            qp.setPen(QColor("#3a6fa8"))
            qp.drawLine(int(x), y0, int(x), y0 + 5)
        qp.setPen(QColor("#bbb"))
        qp.drawLine(pad_l, y1 + 1, pad_l + plot_w, y1 + 1)

    def _draw_primer_marks(self, qp, n, pad_l, pad_r, pad_t, plot_w, plot_h):
        """Outline the primer regions (excluding the GC clamp) and label them
        FP / RP so the primers are not confused with the plain clamp."""
        fp_reg, rp_reg = self.primer_regions()
        pen = QPen(QColor("#3a6fa8"))
        pen.setWidthF(1.1)
        for label, (a, b) in (("FP", fp_reg), ("RP", rp_reg)):
            if b <= a:
                continue
            xa = self._xf(a, pad_l, plot_w, n)
            xb = self._xf(max(b - 1, a), pad_l, plot_w, n) + 1
            qp.setPen(pen)
            qp.drawLine(int(xa), int(pad_t + 2), int(xb), int(pad_t + 2))
            qp.drawLine(int(xa), int(pad_t + plot_h - 2),
                        int(xb), int(pad_t + plot_h - 2))
            qp.setPen(QColor("#3a6fa8"))
            qp.setFont(QFont("Helvetica", 7, QFont.Bold))
            qp.drawText(int(xa + 4), int(pad_t + 3), 26, 10,
                        Qt.AlignLeft, label)

    # -- painting --------------------------------------------------------- #
    def paintEvent(self, event):                        # noqa: N802
        qp = QPainter(self)
        qp.setRenderHint(QPainter.Antialiasing)
        w, h = self.width(), self.height()
        qp.fillRect(0, 0, w, h, QColor("white"))
        pad_l, pad_r, pad_t, pad_b = 46, 12, 30, 26
        plot_w, plot_h = w - pad_l - pad_r, h - pad_t - pad_b

        if not self._ref:
            qp.setPen(QColor("#999"))
            qp.drawText(self.rect(), Qt.AlignCenter, self._status or "no data")
            return

        n = len(self._ref)

        # base-number ruler + primer region markers
        self._draw_ruler(qp, n, pad_l, pad_r, plot_w)

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

        # wt minus mutant melt separation (wt/mut pairs only)
        self._draw_delta_band(qp, n, pad_l, pad_t, plot_w, plot_h)

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

        # primer regions (the GC clamp is not part of the primer)
        self._draw_primer_marks(qp, n, pad_l, pad_r, pad_t, plot_w, plot_h)

        # footer caption
        qp.setFont(QFont("Helvetica", 8))
        qp.setPen(QColor("#999"))
        if self._paired and self._alt:
            cap = ("black = wildtype    dashed red = mutant    "
                   "shaded = melt separation    blue = primers    (5'->3')")
        elif self._alt:
            cap = ("black = reference    dashed red = variant    "
                   "blue = primers    (5'->3')")
        else:
            cap = "reference melt map    blue = primers    (5'->3')"
        qp.drawText(int(pad_l), h - 20, plot_w, 16, Qt.AlignLeft, cap)
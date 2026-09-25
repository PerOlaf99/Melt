"""QPainter-based melting-profile chart for the varmelt GUI.

One canvas overlays the per-base melting maps of all ticked amplicons so
they can be compared under the same conditions (CTCE).  The x-axis is the
*physical fragment*: base 1 is the first base of the DNA that melts, so a
``5'`` GC-clamp oligo occupies bases 1..42 (grey tail) and the amplicon
starts at base 43, while a ``3'`` clamp sits after the amplicon.  Every
dataset is numbered from its own base 1 which is always the left edge of
the plot.

For a wildtype/mutant pair the two GC clamps are the *same* oligo, so with
a ``3'`` clamp the mutant is right-anchored: its clamp is drawn directly
under the wildtype clamp (and the delta band only compares amplicon bases)
so a deletion/insertion can never spill phantoms into the grey clamp tail.

Each item gets its own colour (wildtype solid, mutant dashed) and a legend
shows which line belongs to which amplicon.  There are no primer or
mutation markers.  The view can be zoomed with the wheel or the +/- /
Reset actions, with a drag box, and panned by shift/middle-dragging
(double click resets).
"""
import math

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QSizePolicy, QWidget

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
           "#17becf", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"]

# the mutant (alternative) curve is always drawn dashed in this red, no
# matter how many wildtype items share the chart, so "red dashed" reads
# as "the variant" at a glance.
MUTANT_COLOR = "#d41a1a"


def amp_start(clamp_side, clamp_len):
    """Substrate index of the first amplicon base for a clamped fragment."""
    return (int(clamp_len) if clamp_side == "5'" else 0)


class MeltChart(QWidget):
    """Chart widget; call :meth:`set_datasets` to (re)plot the overlay."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(620, 260)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setAutoFillBackground(True)
        pal = self.palette()
        pal.setColor(self.backgroundRole(), QColor("white"))
        self.setPalette(pal)
        self._datasets = []
        self._status = "tick an amplicon in the list to plot its melt map"
        self._data_x0, self._data_x1 = 0.0, 1.0
        self._data_y0, self._data_y1 = 40.0, 100.0
        self._view_x0, self._view_x1 = 0.0, 1.0
        self._view_y0, self._view_y1 = 40.0, 100.0
        self._legend_lines = 0
        self._delta_count = 0
        self._rubber = None
        self._pan = None
        self._pad_l, self._pad_r = 52, 18
        self._pad_t, self._pad_b = 14, 30

    # --------------------------------------------------------- data input - #
    def set_datasets(self, datasets):
        """Plot a list of datasets; see ``_compute_bounds`` for the fields."""
        self._datasets = list(datasets or [])
        self._compute_bounds()
        self.fit()
        self.update()

    def clear(self, status="no data"):
        self._datasets = []
        self._status = status
        self.update()

    def _compute_bounds(self):
        xs, ys = [], []
        for d in self._datasets:
            ref = d.get("ref_prof") or []
            alt = d.get("alt_prof") or []
            # 3' clamps right-anchor the mutant so its GC clamp overlays
            # the wildtype clamp (see alt_shift); the shorter fragment may
            # not reach x 0.
            shift = 0.0
            if d.get("paired") and alt and d.get("clamp_side") == "3'":
                shift = float(d.get("alt_shift", 0) or 0)
            xs.append((min(0.0, shift),
                       max(len(ref), len(alt) + shift)))
            for t in list(ref) + list(alt):
                if t is not None and not math.isnan(t):
                    ys.append(float(t))
        self._data_x0, self._data_x1 = 0.0, max((b for _, b in xs),
                                                default=100.0)
        if self._data_x1 - self._data_x0 < 1:
            self._data_x1 = self._data_x0 + 1
        y0 = min(ys, default=50.0)
        y1 = max(ys, default=90.0)
        pad = (y1 - y0) * 0.08 or 2.0
        self._data_y0, self._data_y1 = y0 - pad, y1 + pad
        self._view_y0, self._view_y1 = self._data_y0, self._data_y1

    @staticmethod
    def _amp_window(profile_len, clamp_side, clamp_len):
        """(lo, hi) 0-based index range of the amplicon bases inside a
        clamped profile; the GC clamp itself is never part of it."""
        cl = max(0, int(clamp_len))
        if clamp_side == "5'":
            return cl, profile_len
        if clamp_side == "3'":
            return 0, profile_len - cl
        return 0, profile_len

    @staticmethod
    def substrate_x(d, i):
        """x position of substrate base *i*: the physical fragment, so the
        coordinate is simply the base's position (base 1 at x == 0)."""
        return float(i)

    def max_fragment_len(self):
        span = 0
        for d in self._datasets:
            ref = d.get("ref_prof") or []
            span = max(span, len(ref))
        return span

    # --------------------------------------------------------- zoom / pan - #
    def fit(self):
        self._view_x0, self._view_x1 = self._data_x0, self._data_x1
        self._view_y0, self._view_y1 = self._data_y0, self._data_y1
        self.update()

    def reset(self):
        self.fit()

    def zoom_in(self):
        self._zoom_xy(0.5, 0.5, 1.5, 1.5)

    def zoom_out(self):
        self._zoom_xy(0.5, 0.5, 1.0 / 1.5, 1.0 / 1.5)

    def _zoom_xy(self, x_frac, y_frac, xmag, ymag):
        """Divide the visible spans by *mag* around the given fractions;
        the view is kept inside the data ranges."""
        x0, x1 = self._view_x0, self._view_x1
        y0, y1 = self._view_y0, self._view_y1
        dx, dy = self._data_x1 - self._data_x0, self._data_y1 - self._data_y0
        spanx, spany = (x1 - x0) / xmag, (y1 - y0) / ymag
        if 2.0 <= spanx <= 5.0 * dx and 1.5 <= spany <= 5.0 * dy:
            cx = x0 + (x1 - x0) * x_frac
            cy = y0 + (y1 - y0) * y_frac
            self._view_x0 = cx - spanx * x_frac
            self._view_x1 = self._view_x0 + spanx
            self._view_y0 = cy - spany * y_frac
            self._view_y1 = self._view_y0 + spany
        self.update()

    def _clamp_view(self, vx0, vy0):
        dx = self._data_x1 - self._data_x0
        spanx = self._view_x1 - self._view_x0
        spanx = min(spanx, dx)
        vx0 = min(max(vx0, self._data_x0), self._data_x1 - spanx)
        dy = self._data_y1 - self._data_y0
        spany = min(self._view_y1 - self._view_y0, dy)
        vy0 = min(max(vy0, self._data_y0), self._data_y1 - spany)
        self._view_x0, self._view_x1 = vx0, vx0 + spanx
        self._view_y0, self._view_y1 = vy0, vy0 + spany

    # -- geometry ---------------------------------------------------------- #
    def _x(self, v, plot_w):
        span = self._view_x1 - self._view_x0
        return self._pad_l + (v - self._view_x0) / span * plot_w

    def _y(self, t, plot_h):
        span = self._view_y1 - self._view_y0
        return self._pad_t + (1.0 - (t - self._view_y0) / span) * plot_h

    # -- drawing helpers --------------------------------------------------- #
    def _draw_polyline(self, qp, points, color, width, dash=False):
        pen = QPen(QColor(color))
        pen.setWidthF(width)
        if dash:
            pen.setStyle(Qt.DashLine)
        qp.setPen(pen)
        seg = []
        for x, y in points:
            if y is None or math.isnan(y):
                if len(seg) >= 2:
                    qp.drawPolyline(QPolygonF(seg))
                seg = []
                continue
            seg.append(QPointF(x, y))
        if len(seg) >= 2:
            qp.drawPolyline(QPolygonF(seg))

    def _draw_delta_band(self, qp, d, plot_w, plot_h):
        """Fill between the ref and alt curve of one pair, in its colour."""
        ref = d.get("ref_prof") or []
        alt = d.get("alt_prof")
        if not ref or not alt or not d.get("paired"):
            return
        var = d.get("mark_idx")
        var = 0 if var is None else int(var)
        indel = int(d.get("indel", 0))
        n_alt = len(alt)
        # only amplicon-vs-amplicon bases are compared: the GC clamp is the
        # same oligo on both, so a deletion/insertion must never show up as
        # phantom differences inside the grey clamp tail
        cl = max(0, int(d.get("clamp_len", 0)))
        r_lo, r_hi = self._amp_window(len(ref), d.get("clamp_side"), cl)
        a_lo, a_hi = self._amp_window(n_alt, d.get("clamp_side"), cl)
        top, bot = [], []

        def flush():
            if len(top) >= 2 and len(bot) >= 2:
                col = QColor(d["color"])
                col.setAlpha(60)
                qp.setPen(Qt.NoPen)
                qp.setBrush(col)
                qp.drawPolygon(QPolygonF(
                    [QPointF(x, y) for x, y in top]
                    + [QPointF(x, y) for x, y in reversed(bot)]))
                self._delta_count += 1

        for i, rt in enumerate(ref):
            j = i if i < var else i - indel
            if not (r_lo <= i < r_hi) \
                    or rt is None or math.isnan(rt) \
                    or not (a_lo <= j < a_hi) or alt[j] is None \
                    or math.isnan(alt[j]):
                flush()
                top, bot = [], []
                continue
            x = self._x(self.substrate_x(d, i), plot_w)
            yr = self._y(rt, plot_h)
            ya = self._y(alt[j], plot_h)
            if yr <= ya:
                top.append((x, yr))
                bot.append((x, ya))
            else:
                top.append((x, ya))
                bot.append((x, yr))
        flush()

    def _draw_legend(self, qp, plot_w, plot_h):
        entries = []
        for d in self._datasets:
            name = d.get("name", "")
            if d.get("dip"):
                name += "   \u26a0  valley"
            entries.append((d["color"], True, name, None))
            if d.get("paired") and d.get("alt_prof"):
                area = float(d.get("delta_area", 0.0))
                entries.append((MUTANT_COLOR, False,
                                f"mutant   \u0394area {area:.0f} \u00b0C\u00b7bp",
                                d.get("name")))
        if not entries:
            self._legend_lines = 0
            self._legend_texts = []
            return
        self._legend_texts = [t for _c, _sol, t, _d in entries]
        qp.setFont(QFont("Helvetica", 8))
        fm = qp.fontMetrics()
        row_h = fm.height() + 4
        padx, pady = 10, 8
        x0 = self._pad_l + plot_w - padx
        y = pady
        # measure the widest entry text
        wmax = 0
        for _c, _sol, text, _d in entries:
            wmax = max(wmax, fm.horizontalAdvance(text or ""))
        bw = padx * 2 + 18 + wmax
        bh = pady * 2 + row_h * len(entries)
        box = QRectF(x0 - bw, y, bw, bh)
        qp.setBrush(QColor(255, 255, 255, 215))
        qp.setPen(QPen(QColor("#ccc")))
        qp.drawRect(box)
        y += pady + fm.ascent() - 1
        sw = 12
        for color, solid, text, _d in entries:
            sx = x0 - bw + padx
            qp.setPen(Qt.NoPen)
            qp.setBrush(Qt.NoBrush)
            pen = QPen(QColor(color))
            pen.setWidthF(1.6)
            if not solid:
                pen.setStyle(Qt.DashLine)
            qp.setPen(pen)
            qp.drawLine(int(sx), int(y - fm.ascent() / 2),
                        int(sx + sw), int(y - fm.ascent() / 2))
            qp.setPen(QColor("#333"))
            qp.drawText(int(sx + sw + 5), int(y - fm.ascent()),
                        fm.horizontalAdvance(text) + 4, fm.height(),
                        Qt.AlignLeft | Qt.AlignVCenter, text)
            y += row_h
        self._legend_lines = len(entries)

    # -- painting --------------------------------------------------------- #
    def paintEvent(self, event):                        # noqa: N802
        qp = QPainter(self)
        qp.setRenderHint(QPainter.Antialiasing)
        w, h = self.width(), self.height()
        qp.fillRect(0, 0, w, h, QColor("white"))
        pad_l, pad_r, pad_t, pad_b = (self._pad_l, self._pad_r,
                                      self._pad_t, self._pad_b)
        plot_w, plot_h = w - pad_l - pad_r, h - pad_t - pad_b

        if not self._datasets:
            qp.setPen(QColor("#999"))
            qp.drawText(self.rect(), Qt.AlignCenter, self._status or "no data")
            return

        self._legend_lines = 0
        self._delta_count = 0
        span_x = self._data_x1 - self._data_x0

        # clip every band and curve to the plot box so zoomed/panned views
        # never let a line cross the axes into the margins
        qp.save()
        qp.setClipRect(QRectF(pad_l, pad_t, plot_w, plot_h))

        # physical fragment band behind each dataset's amplicon (base 1..N),
        # with the GC clamp drawn as a grey tail at the end that carries it
        for d in self._datasets:
            ref = d.get("ref_prof") or []
            cl = max(0, int(d.get("clamp_len", 0)))
            side = d.get("clamp_side")
            amp = max(len(ref) - cl, 0)
            a = amp_start(side, cl)
            if amp > 0:
                qp.setPen(Qt.NoPen)
                qp.setBrush(QColor("#f7fafd"))
                x0 = self._x(a, plot_w)
                x1 = self._x(a + amp, plot_w)
                qp.drawRect(int(x0), int(pad_t), int(max(x1 - x0, 2)),
                            int(plot_h))
            if cl and side in ("5'", "3'"):
                ca = 0 if side == "5'" else amp
                qp.setPen(Qt.NoPen)
                qp.setBrush(QColor("#efeceb"))
                x0 = self._x(ca, plot_w)
                x1 = self._x(ca + cl, plot_w)
                qp.drawRect(int(x0), int(pad_t), int(max(x1 - x0, 2)),
                            int(plot_h))
                qp.setPen(QColor("#b6a9a4"))
                qp.setFont(QFont("Helvetica", 7, QFont.Bold))
                qp.drawText(int((x0 + x1) / 2 - 46), int(pad_t + 5), 92, 12,
                            Qt.AlignHCenter, f"GC clamp {side}")

        # per-pair wt/mut separation band in that item's colour
        for d in self._datasets:
            self._draw_delta_band(qp, d, plot_w, plot_h)

        # curves: solid reference / dashed mutant, one colour per item
        qp.setBrush(Qt.NoBrush)
        for d in self._datasets:
            ref = d.get("ref_prof") or []
            color = d["color"]
            pts = [(self._x(self.substrate_x(d, i), plot_w),
                    self._y(t, plot_h)) for i, t in enumerate(ref)]
            self._draw_polyline(qp, pts, color, 1.6)
            alt = d.get("alt_prof") or []
            if alt and d.get("paired"):
                # right-anchored for 3' clamps: keep the identical GC clamp
                # glued to the wildtype clamp, so an indel never appears to
                # extend into the grey clamp tail
                alt_shift = float(d.get("alt_shift", 0) or 0)
                pts = [(self._x(self.substrate_x(d, i) + alt_shift, plot_w),
                        self._y(t, plot_h)) for i, t in enumerate(alt)]
                self._draw_polyline(qp, pts, MUTANT_COLOR, 1.2, dash=True)
        qp.restore()

        # axes: y grid + labels, x grid labelled with physical base numbers
        qp.setFont(QFont("Helvetica", 8))
        step = max(int(math.ceil((self._data_y1 - self._data_y0) / 12)), 5)
        for t in range(int(self._data_y0 // step) * step,
                       int(self._data_y1) + 1, step):
            if t < self._data_y0 or t > self._data_y1:
                continue
            y = self._y(t, plot_h)
            qp.setPen(QColor("#e4e4e4"))
            qp.drawLine(int(pad_l), int(y), int(w - pad_r), int(y))
            qp.setPen(QColor("#888"))
            qp.drawText(int(pad_l - 46), int(y - 7), 42, 14,
                        Qt.AlignRight | Qt.AlignVCenter, f"{t:g} °C")

        qp.setPen(QColor("#a8a8a8"))
        qp.drawLine(int(pad_l), int(pad_t), int(pad_l),
                    int(pad_t + plot_h))
        qp.drawLine(int(pad_l), int(pad_t + plot_h),
                    int(w - pad_r), int(pad_t + plot_h))

        qp.setFont(QFont("Helvetica", 7))
        tick = 20 if span_x > 180 else 10 if span_x > 90 else 5
        for b in range(0, int(span_x) + 1, tick):
            x = self._x(b, plot_w)
            if not (pad_l - 1 <= x <= w - pad_r + 1):
                continue
            qp.setPen(QColor("#cfcfcf"))
            qp.drawLine(int(x), int(pad_t + plot_h),
                        int(x), int(pad_t + plot_h + 3))
            qp.setPen(QColor("#777"))
            qp.drawText(int(x - 14), int(pad_t + plot_h + 2), 28, 12,
                        Qt.AlignHCenter, str(b + 1))
        qp.setPen(QColor("#999"))
        qp.drawText(int(w - pad_r - 240), h - 13, 240, 12,
                    Qt.AlignRight,
                    "drag box zoom · shift/mid-drag pan · wheel/± zoom · "
                    "double-click reset")

        self._draw_legend(qp, plot_w, plot_h)
        if self._rubber is not None:
            self._draw_rubber(qp, plot_w, plot_h)

    def _draw_rubber(self, qp, plot_w, plot_h):
        r = self._rubber
        w, h = self.width(), self.height()
        x = max(r.left(), self._pad_l)
        y = max(r.top(), self._pad_t)
        x2 = min(r.right(), w - self._pad_r)
        y2 = min(r.bottom(), self._pad_t + plot_h)
        if x2 - x < 1 or y2 - y < 1:
            return
        qp.setBrush(QColor(90, 160, 220, 55))
        qp.setPen(QPen(QColor(50, 120, 190), 1.2))
        qp.drawRect(int(x), int(y), int(x2 - x), int(y2 - y))

    # ---------------------------------------------------------- interactions - #
    def wheelEvent(self, event):                        # noqa: N802
        pos = event.position()
        plot_w = max(self.width() - self._pad_l - self._pad_r, 1)
        plot_h = max(self.height() - self._pad_t - self._pad_b, 1)
        x_frac = min(max((pos.x() - self._pad_l) / plot_w, 0.0), 1.0)
        y_frac = min(max((pos.y() - self._pad_t) / plot_h, 0.0), 1.0)
        delta = event.angleDelta().y()
        if delta:
            mag = 1.35 if delta > 0 else 1.0 / 1.35
            self._zoom_xy(x_frac, y_frac, mag, mag)
        event.accept()

    def mousePressEvent(self, event):                   # noqa: N802
        if event.button() == Qt.LeftButton:
            if event.modifiers() & Qt.ShiftModifier:
                # Shift+drag pans (like the middle button); plain left-drag
                # rubber-bands a zoom rectangle.
                self._rubber = None
                self._pan = (event.position().x(), event.position().y(),
                             self._view_x0, self._view_x1,
                             self._view_y0, self._view_y1)
            else:
                self._pan = None
                px, py = event.position().x(), event.position().y()
                self._rubber = QRectF(px, py, 0.0, 0.0)
            event.accept()
        elif event.button() == Qt.MiddleButton:
            self._rubber = None
            self._pan = (event.position().x(), event.position().y(),
                         self._view_x0, self._view_x1,
                         self._view_y0, self._view_y1)
            event.accept()

    def mouseMoveEvent(self, event):                    # noqa: N802
        if self._pan is not None:
            x_0, y_0, vx0, vx1, vy0, vy1 = self._pan
            plot_w = max(self.width() - self._pad_l - self._pad_r, 1)
            plot_h = max(self.height() - self._pad_t - self._pad_b, 1)
            spanx, spany = vx1 - vx0, vy1 - vy0
            shift_x = spanx * (event.position().x() - x_0) / plot_w
            shift_y = spany * (event.position().y() - y_0) / plot_h
            self._clamp_view(vx0 - shift_x, vy0 + shift_y)
        elif self._rubber is not None:
            r = self._rubber
            r.setRight(event.position().x())
            r.setBottom(event.position().y())
            self._rubber = r.normalized()
            self.update()

    def mouseReleaseEvent(self, event):                 # noqa: N802
        if event.button() == Qt.LeftButton \
                and self._rubber is not None:
            rect = self._rubber
            self._rubber = None
            self.update()
            if rect.width() >= 6 and rect.height() >= 6:
                self._zoom_rect(rect)
        elif self._pan is not None:
            self._pan = None

    def keyPressEvent(self, event):                     # noqa: N802
        if event.key() == Qt.Key_Escape and self._rubber is not None:
            self._rubber = None
            self.update()
            event.accept()
        else:
            super().keyPressEvent(event)

    def mouseDoubleClickEvent(self, event):             # noqa: N802
        self._rubber = None
        self.fit()
        event.accept()

    def _zoom_rect(self, rect):
        """Zoom the plot so the on-screen rectangle *rect* (pixel coords,
        plot-area coordinate space) fills the view.  Tiny rectangles are
        rejected so an accidental click never collapses the axis."""
        plot_w = max(self.width() - self._pad_l - self._pad_r, 1)
        plot_h = max(self.height() - self._pad_t - self._pad_b, 1)
        pad_l, pad_t = self._pad_l, self._pad_t
        x = min(max(rect.left(), pad_l), pad_l + plot_w)
        x2 = min(max(rect.right(), pad_l), pad_l + plot_w)
        y = min(max(rect.top(), pad_t), pad_t + plot_h)
        y2 = min(max(rect.bottom(), pad_t), pad_t + plot_h)
        if x2 - x < 6 or y2 - y < 6:
            return
        spanx = self._view_x1 - self._view_x0
        spany = self._view_y1 - self._view_y0
        # pixel y grows downwards while the map's hot end sits at the top
        # (small y): the box top is the high-Tm edge, the bottom the low end.
        vx = self._view_x0 + (x - pad_l) / plot_w * spanx
        t_top = self._view_y0 + (plot_h - (y - pad_t)) / plot_h * spany
        t_bot = self._view_y0 + (plot_h - (y2 - pad_t)) / plot_h * spany
        new_spanx = spanx * (x2 - x) / plot_w
        new_spany = spany * (y2 - y) / plot_h
        if new_spanx < 2.0 or new_spany < 1.5:
            return
        self._view_x0, self._view_x1 = vx, vx + new_spanx
        self._view_y0, self._view_y1 = t_bot, t_top
        self._clamp_view(self._view_x0, self._view_y0)
        self.update()
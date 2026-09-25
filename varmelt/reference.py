"""Reference DNA-melting model faithful to the original tool.

This is a port of the Perl package `DNAMelting` (the implementation used
by the *Genomic HyperBrowser* "Variant melting profiles" tool / MeltPrimer,
by Tostesen, Nakken et al.).  Reference: E. Tøstesen, F. Liu, T.-K. Jenssen,
E. Hovig, *Biopolymers* 70:364-376 (2003).

Thermodynamics  : Blake & Delcourt (1998) nearest-neighbour parameters as
                  reparameterised by Blossey & Carlon (2003): each stack has
                  its own melting temperature ``Tij`` that depends on the
                  monovalent salt concentration.

Partition model: Yeramian & Tøstesen forward/backward recursions with a
                  multiexponential approximation of the exact loop-entropy
                  factor ``Omega(2k) = sigma * (2k+d)^(-alpha)``
                  (``sigma = 1.26e-4``, ``alpha = 2.15``), giving O(N) time
                  per temperature.

Profiles:
  * ``calc_tm_profile``  : per-base local melting temperature (deg C), the
                           temperature at which the closed probability of
                           each base crosses ``helicity`` (default 0.5).
  * ``melt_curve``       : fraction of closed (paired) bases vs temperature,
                           i.e. the full melting curve of the sequence.
"""

import math

# --- Blake-Delcourt 1998 parameters (reparameterised, Blossey-Carlon 03) --
# Keyed by the canonical dinucleotide names used in the Perl module.
T0ij = {
    "AA=TT": 89.08, "AT": 89.38, "TA": 79.47, "CA=TG": 89.71,
    "AC=GT": 121.17, "AG=CT": 98.49, "GA=TC": 105.09, "CG": 105.28,
    "GC": 143.73, "CC=GG": 118.49,
}

dHij = {
    "AC=GT": 10.51, "TA": 7.81, "AA=TT": 8.45, "GC": 11.91,
    "CC=GG": 10.34, "GA=TC": 9.47, "AT": 8.50, "CG": 9.53,
    "AG=CT": 9.10, "CA=TG": 8.51,
}

dTij_dlog10Na = {
    "AC=GT": 13.71, "TA": 22.08, "AA=TT": 19.78, "GC": 10.21,
    "CC=GG": 14.18, "GA=TC": 16.94, "AT": 19.19, "CG": 16.01,
    "AG=CT": 17.30, "CA=TG": 19.35,
}

nn2ten = {
    "AA": "AA=TT", "TT": "AA=TT", "AT": "AT", "TA": "TA",
    "CA": "CA=TG", "TG": "CA=TG", "GT": "AC=GT", "AC": "AC=GT",
    "CT": "AG=CT", "AG": "AG=CT", "GA": "GA=TC", "TC": "GA=TC",
    "CG": "CG", "GC": "GC", "GG": "CC=GG", "CC": "CC=GG",
}

R = 1.987      # gas constant, cal/(mol K)
d = 1          # loop-diameter correction
SIGMA = 1.26e-4
ALPHA = 2.15
BETA = 1.0

_LOG10E = 1.0 / math.log(10.0)


# ---------------------------------------------------------------------------
# Loop-entropy factor: multiexponential approximation (Yeramian-Tostesen)
# ---------------------------------------------------------------------------
def _muex(n: int):
    """Return tuples ``(C1, C2)`` for the multiexponential loop factor.

    Sets ``Omega(2k) = sum_m C1[m] * C2[m]**(k-2)`` for ``k = 2..n``.
    """
    I = int(1 + math.log(2 * n))
    B = [0.0] * (I + 1)
    A = [0.0] * (I + 1)
    for m in range(I + 1):
        B[m] = math.exp(m - I)
        A[m] = math.exp(1 - ALPHA * (I - m))
        for mm in range(m):
            A[m] -= A[mm] * math.exp(1 - math.exp(mm - m))
    k = sum(A[m] * math.exp(-B[m] * math.exp(3)) for m in range(I + 1))
    for m in range(I + 1):
        A[m] /= k * math.exp(ALPHA * 3)

    C1 = [SIGMA * A[m] * math.exp(-d * B[m]) * math.exp(-4.0 * B[m])
          for m in range(I + 1)]
    C2 = [math.exp(-2.0 * B[m]) for m in range(I + 1)]
    return C1, C2


# ---------------------------------------------------------------------------
# Statistical weights (Blossey-Carlon / Blake-Delcourt)
# ---------------------------------------------------------------------------
def _stack_keys(seq: str):
    """Dinucleotide keys (canonical) of the sequence, indexed 2..N."""
    keys = {}
    for i in range(1, len(seq)):
        keys[i + 1] = nn2ten[seq[i - 1:i + 1]]
    return keys


def _s_nn(keys: dict, T_K: float, Na: float):
    """Forward (``s_LR[i]``) and reversed (``s_RL``) stack weights.

    ``sij = exp( -1000*dHij/(R*Tij) * (1 - Tij/T) )`` with
    ``Tij = T0ij + 273.15 + log10(Na)*dTij`` (T0ij given in deg C).
    """
    log10Na = math.log(Na) * _LOG10E
    Tij = {key: T0ij[key] + 273.15 + log10Na * dTij_dlog10Na[key]
           for key in keys.values()}
    N = max(keys)          # sequence length (index of last stack + 1)
    s_LR = [0.0] * (N + 2)
    for i, key in keys.items():
        T = Tij[key]
        s_LR[i] = math.exp((-1000.0 * dHij[key] / (R * T)) * (1.0 - T / T_K))
    s_RL = [0.0] * (N + 2)
    for i in range(2, N + 1):
        s_RL[N + 2 - i] = s_LR[i]
    return s_LR, s_RL


# ---------------------------------------------------------------------------
# Partition function (O(N), Yeramian-Tostesen recursions)
# ---------------------------------------------------------------------------
def _partition(s_LR, s_RL, C1, C2):
    """Compute LR/RL recursion arrays and total partition function.

    Returns ``(V_LR, U1_LR, U2_LR, Z_LR, V_RL, U1_RL, U2_RL)`` all indexed
    by 1..N+1 in the Perl convention.
    """
    N = len(s_LR) - 2
    M = len(C1) - 1

    norm_factor = 1e-30
    bignumber = 1e60

    norm_times_LR = {}
    norm_times_RL = {}
    norm_times_RL = {}

    # ---- Left-Right (forward) ----
    V_LR = [0.0] * (N + 3)
    U1_LR = [0.0] * (N + 3)
    U2_LR = [0.0] * (N + 3)
    V_LR[1] = 1.0
    V_LR[2] = BETA * 0.0                      # beta * s_solo_LR[1] = 0
    U1_LR[1] = 1.0
    U1_LR[2] = BETA
    U2_LR[2] = BETA * 1.0 * s_LR[2]           # beta * s_end_LR[1] * s_NN_LR[2]
    V_LR[3] = U2_LR[2] * 1.0 + U1_LR[2] * 0.0  # s_end=1, s_solo=0
    W = [V_LR[2] * C1[m] for m in range(M + 1)]
    Z_LR = V_LR[1] + V_LR[2] + V_LR[3]
    norm_times_RL = {}

    for i in range(3, N + 1):
        U1_LR[i] = BETA * V_LR[1] + sum(W)
        U2_LR[i] = s_LR[i] * (U1_LR[i - 1] * 1.0 + U2_LR[i - 1])
        V_LR[i + 1] = U2_LR[i] * 1.0 + U1_LR[i] * 0.0
        for m in range(M + 1):
            W[m] = V_LR[i] * C1[m] + W[m] * C2[m]
        Z_LR += V_LR[i + 1]
        if Z_LR > bignumber:
            Z_LR *= norm_factor
            U1_LR[i] *= norm_factor
            U2_LR[i] *= norm_factor
            V_LR[i + 1] *= norm_factor
            V_LR[1] *= norm_factor
            for m in range(M + 1):
                W[m] *= norm_factor
            norm_times_RL[N + 2 - i] = 1

    # ---- Right-Left (backward) ----
    V_RL = [0.0] * (N + 3)
    U1_RL = [0.0] * (N + 3)
    U2_RL = [0.0] * (N + 3)
    V_RL[1] = 1.0
    V_RL[2] = BETA * 0.0
    U1_RL[1] = 1.0
    U1_RL[2] = BETA
    U2_RL[2] = BETA * 1.0 * s_RL[2]
    V_RL[3] = U2_RL[2] * 1.0 + U1_RL[2] * 0.0
    W = [V_RL[2] * C1[m] for m in range(M + 1)]
    Z_RL = V_RL[1] + V_RL[2] + V_RL[3]

    if norm_times_RL.get(2):
        Z_RL *= norm_factor
        U1_RL[2] *= norm_factor
        U2_RL[2] *= norm_factor
        V_RL[3] *= norm_factor
        V_RL[1] *= norm_factor
        W = [w * norm_factor for w in W]

    for i in range(3, N + 1):
        U1_RL[i] = BETA * V_RL[1] + sum(W)
        U2_RL[i] = s_RL[i] * (U1_RL[i - 1] * 1.0 + U2_RL[i - 1])
        V_RL[i + 1] = U2_RL[i] * 1.0 + U1_RL[i] * 0.0
        for m in range(M + 1):
            W[m] = V_RL[i] * C1[m] + W[m] * C2[m]
        Z_RL += V_RL[i + 1]
        if norm_times_RL.get(i):
            Z_RL *= norm_factor
            U1_RL[i] *= norm_factor
            U2_RL[i] *= norm_factor
            V_RL[i + 1] *= norm_factor
            V_RL[1] *= norm_factor
            W = [w * norm_factor for w in W]

    return (V_LR, U1_LR, U2_LR, Z_LR, V_RL, U1_RL, U2_RL, Z_RL)


def _p_closed(i, V_LR, U1_LR, U2_LR, Z_LR, V_RL, U1_RL, U2_RL, Z_RL, N):
    """Closed probability of base ``i`` (1-based), Perl p_closed()."""
    if i == 0:
        return 0.0
    if i == 1:
        return V_RL[N + 1] / Z_RL
    if i == N:
        return V_LR[N + 1] / Z_LR
    return ((U1_LR[i] * U2_RL[N + 1 - i]
             + U2_LR[i] * U1_RL[N + 1 - i]
             + U2_LR[i] * U2_RL[N + 1 - i]) / (BETA * Z_LR))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def calc_tm_profile(seq: str, Na: float = 0.013, helicity: float = 0.5,
                    tmin: int = 40, tmax: int = 120) -> list:
    """Per-base local melting temperature (deg C), one value per base.

    Port of ``Meltingprofiles::calc_temp_profile``: scans temperature from
    ``tmin`` to ``tmax`` deg C (1 K steps) and records, for every base, the
    temperature at which the closed probability crosses ``helicity``,
    interpolated in logit/1-T space (sigmoidal) or linearly for p~{0,1}.
    """
    seq = seq.upper()
    N = len(seq)
    if N < 3:
        return [float('nan')] * N

    keys = _stack_keys(seq)
    C1, C2 = _muex(N)

    T = 59.0 + 273.0
    s_LR, s_RL = _s_nn(keys, T, Na)
    arrays = _partition(s_LR, s_RL, C1, C2)
    p1 = [_p_closed(i, *arrays, N) for i in range(1, N + 1)]

    TL = [float('nan')] * (N + 1)      # 1-based; TL[0] unused ('' in Perl)
    pL = helicity
    FL = math.log(pL / (1.0 - pL))

    for T in range(tmin + 273, tmax + 274):     # deg C -> K, inclusive
        TI = T / 1.0
        s_LR, s_RL = _s_nn(keys, TI, Na)
        arrays = _partition(s_LR, s_RL, C1, C2)
        for i in range(1, N + 1):
            p2 = _p_closed(i, *arrays, N)
            if (p1[i - 1] - pL) * (p2 - pL) <= 0:
                p_prev = p1[i - 1]
                if abs(p_prev - 0.5) < 0.5 and abs(p2 - 0.5) < 0.5:
                    F1 = math.log(p_prev / (1.0 - p_prev))
                    F2 = math.log(p2 / (1.0 - p2))
                    V1 = 1.0 / (TI - 1.0)
                    V2 = 1.0 / TI
                    VL = (V1 * ((FL - F2) / (F1 - F2))
                          + V2 * ((FL - F1) / (F2 - F1)))
                    TL[i] = 1.0 / VL - 273.0
                else:
                    TL[i] = (-273.0 + (TI - 1.0) * ((pL - p2) / (p_prev - p2))
                             + TI * ((pL - p_prev) / (p2 - p_prev)))
            p1[i - 1] = p2

    return [TL[i] for i in range(1, N + 1)]


def melt_curve(seq: str, Na: float = 0.013, temps=None) -> list:
    """Fraction of closed (base-paired) bases vs temperature (deg C).

    Returns a list of ``(T_degC, p_closed)`` points; ``melting`` is the
    temperature at which half the bases are paired.
    """
    seq = seq.upper()
    N = len(seq)
    if N < 3:
        return []
    keys = _stack_keys(seq)
    C1, C2 = _muex(N)
    if temps is None:
        temps = [t for t in range(20, 101, 1)]

    out = []
    for T_deg in temps:
        T = T_deg + 273.15
        s_LR, s_RL = _s_nn(keys, T, Na)
        arrays = _partition(s_LR, s_RL, C1, C2)
        p = sum(_p_closed(i, *arrays, N) for i in range(1, N + 1)) / N
        out.append((T_deg, p))
    return out


def melting_temp(seq: str, Na: float = 0.013, helicity: float = 0.5,
                 tmin: int = 40, tmax: int = 120) -> float:
    """Whole-sequence melting temperature (deg C) at the given helicity."""
    prof = calc_tm_profile(seq, Na=Na, helicity=helicity,
                           tmin=tmin, tmax=tmax)
    vals = [v for v in prof if not math.isnan(v)]
    return sum(vals) / len(vals) if vals else float('nan')
"""Light-curve helpers shared across the package.

These are deliberately minimal: a CSV loader with faint-outlier clipping, a
cubic-polynomial detrender, and a K-harmonic Fourier peak-to-peak amplitude.
Input format: CSV with columns ``time`` (BTJD days), ``mag``, and optionally
``err``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# TESS spacecraft orbit (hours): momentum-dump systematics recur at ORBIT_H / n.
ORBIT_H = 328.8


def load_clean(path: str):
    """Load a (time, mag, err) CSV; drop non-finite rows and iteratively clip
    BRIGHT outliers (m < median - 4*MAD, e.g. background-star spikes).
    Returns (t, mag, err, n_clipped), time-sorted."""
    df = pd.read_csv(path)
    t = pd.to_numeric(df["time"], errors="coerce").values
    m = pd.to_numeric(df["mag"], errors="coerce").values
    e = pd.to_numeric(df["err"], errors="coerce").values if "err" in df.columns else np.full(len(t), 0.1)
    g = np.isfinite(t) & np.isfinite(m) & np.isfinite(e) & (e < 1.0)
    if g.sum() < 20:
        g = np.isfinite(t) & np.isfinite(m)
    t, m, e = t[g], m[g], e[g]
    n_clip = 0
    for _ in range(3):
        med = np.median(m)
        mad = 1.4826 * np.median(np.abs(m - med))
        if mad <= 0:
            break
        bad = m < med - 4 * mad
        if not bad.any():
            break
        n_clip += int(bad.sum())
        t, m, e = t[~bad], m[~bad], e[~bad]
    o = np.argsort(t)
    return t[o], m[o], e[o], n_clip


def detrend(t: np.ndarray, m: np.ndarray, deg: int = 3) -> np.ndarray:
    """Remove a degree-``deg`` polynomial trend; returns zero-mean residuals."""
    x = (t - t.mean()) / max(t.max() - t.min(), 1e-9)
    return m - np.polyval(np.polyfit(x, m, deg), x)


def fourier_amp(t, y, P_rot_h, K=2):
    """Peak-to-peak of a K-harmonic Fourier fit at P_rot_h (t in days).
    Scatter-robust. Returns (peak_to_peak, A1, A2)."""
    w = 2 * np.pi / (P_rot_h / 24.0)
    cols = [np.ones_like(t)]
    for k in range(1, K + 1):
        cols += [np.cos(k * w * t), np.sin(k * w * t)]
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), y, rcond=None)
    ph = np.linspace(0, 1, 400)
    tt = ph * (P_rot_h / 24.0)
    c2 = [np.ones_like(tt)]
    for k in range(1, K + 1):
        c2 += [np.cos(k * w * tt), np.sin(k * w * tt)]
    mm = np.column_stack(c2) @ coef
    A1 = float(np.hypot(coef[1], coef[2]))
    A2 = float(np.hypot(coef[3], coef[4])) if K >= 2 else 0.0
    return float(mm.max() - mm.min()), A1, A2

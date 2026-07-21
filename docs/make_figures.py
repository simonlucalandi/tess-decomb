#!/usr/bin/env python3
"""
Regenerate the explainer figures embedded in the top-level README.

Every panel is drawn from real data committed to this repository:
  fig1_comb.png      analytic comb (328.8/n h) + three ground-truth cases
  fig2_plane.png     validation/survey_verdicts.csv   (221 real verdicts)
  fig3_injection.png validation/v2_injections.csv     (17,372 injections)
  fig4_validity.png  validation/v6_validity.csv       (validity domain)

Usage:  python3 docs/make_figures.py
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
VAL = ROOT / "validation"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)

DUMP_H = 328.8            # momentum-dump comb fundamental
SURVIVE_DROP = 15.0       # <= this: survived
KILL_DROP = 50.0          # >= this AND r2 >= KILL_R2: killed
KILL_R2 = 0.30

INK = "#101720"
DIM = "#4A5A67"
GRID = "#D5DEE4"
SIG = "#1F7A63"           # survived / astrophysical
INST = "#A8382A"          # killed / instrumental
AMB = "#A2761A"           # inconclusive

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": DIM, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": DIM, "ytick.color": DIM,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  {path.relative_to(ROOT)}")


def fig_comb():
    """The alias comb, and why proximity to a tooth is not a verdict."""
    fig, ax = plt.subplots(figsize=(8.6, 2.5))
    for n in range(1, 13):
        p = DUMP_H / n
        if not 20 <= p <= 400:
            continue
        ax.axvline(p, color=GRID, lw=6 if n <= 6 else 4, zorder=1)
        if n <= 6:
            ax.text(p, 1.06, f"n={n}", ha="center", va="bottom",
                    fontsize=7.5, color=DIM, family="monospace")
    # ground-truth cases from the K-sweep validation (see validation/V1_RESULTS.md)
    cases = [(108.08, "(14720) 108.08 h\nreal rotation", SIG, 0.72),
             (85.00, "(1124) 85.0 h\ncomb alias", INST, 0.42),
             (25.88, "(14274) 25.88 h\nunresolved", AMB, 0.72)]
    for p, label, col, y in cases:
        ax.plot([p, p], [0, y], color=col, lw=1.8, zorder=3)
        ax.plot([p], [y], "o", color=col, ms=6, zorder=4)
        ha = "right" if p > 200 else ("left" if p < 40 else "center")
        ax.text(p, y + 0.06, label, ha=ha, va="bottom", fontsize=7.8, color=col,
                linespacing=1.35)
    ax.set_xscale("log")
    ax.set_xlim(20, 400)
    ax.set_ylim(0, 1.30)
    ax.set_yticks([])
    ax.set_xticks([20, 30, 50, 80, 120, 200, 330])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("period (hours, log scale)")
    ax.spines["left"].set_visible(False)
    ax.set_title("TESS momentum-dump alias comb at 328.8/n h", loc="left",
                 fontsize=9.5, color=INK, pad=14)
    save(fig, "fig1_comb.png")


def fig_plane():
    """Where the decision rule actually places 221 real survey verdicts."""
    rows = list(csv.DictReader(open(VAL / "survey_verdicts.csv")))
    by = defaultdict(list)
    for r in rows:
        by[r["verdict"]].append((float(r["drop_pct"]), float(r["r2"]), r["number"],
                                 float(r["period_h"])))

    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ax.axvspan(-12, SURVIVE_DROP, color=SIG, alpha=0.09, zorder=0)
    ax.add_patch(plt.Rectangle((KILL_DROP, KILL_R2), 100, 10,
                               color=INST, alpha=0.09, zorder=0))
    ax.axvline(SURVIVE_DROP, color=SIG, ls="--", lw=1.2)
    ax.axvline(KILL_DROP, color=INST, ls="--", lw=1.2)
    ax.plot([KILL_DROP, 75], [KILL_R2, KILL_R2], color=INST, ls="--", lw=1.2)

    style = {"survived": (SIG, 22, 0.45), "inconclusive": (AMB, 30, 0.85),
             "KILLED": (INST, 130, 1.0)}
    for verdict, (col, size, alpha) in style.items():
        pts = by.get(verdict, [])
        if not pts:
            continue
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=size, c=col,
                   alpha=alpha, linewidths=0, zorder=3,
                   label=f"{verdict}  n={len(pts)}")
    for d, r2, num, per in by.get("KILLED", []):
        ax.annotate(f"({num})  {per:g} h", (d, r2), xytext=(-14, 0),
                    textcoords="offset points", ha="right", va="center",
                    fontsize=8.5, color=INST, family="monospace")
    # the instructive near-miss: biggest drop in the whole run, but the
    # systematics model explained essentially nothing, so the R^2 gate spares it
    for d, r2, num, per in by.get("inconclusive", []):
        if d > KILL_DROP:
            ax.annotate(f"({num}) {per:g} h\ndrop {d:g}% but $R^2$={r2:g}\n"
                        "-> spared by the $R^2$ gate",
                        (d, r2), xytext=(-6, 30), textcoords="offset points",
                        ha="right", va="bottom", fontsize=7.6, color=DIM,
                        linespacing=1.4,
                        arrowprops=dict(arrowstyle="-", color=DIM, lw=0.8))

    ax.set_xlim(-12, 78)
    ax.set_ylim(0, 0.92)
    ax.set_xlabel("Lomb-Scargle power drop after projection  (%)")
    ax.set_ylabel("systematics-model $R^2$")
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.5, loc="lower center",
              bbox_to_anchor=(0.5, 1.02), ncol=3)
    ax.set_title("221 real verdicts from the survey runs", loc="left",
                 fontsize=9.5, color=INK, pad=34)
    save(fig, "fig2_plane.png")


def fig_injection():
    """False-kill rate by injected amplitude and placement relative to a tooth."""
    tally = defaultdict(lambda: [0, 0])           # (kind, amp) -> [killed, total]
    for r in csv.DictReader(open(VAL / "v2_injections.csv")):
        key = (r["kind"], float(r["amp_in"]))
        tally[key][1] += 1
        if r["verdict"] == "killed":
            tally[key][0] += 1
    kinds = [("tooth", "on a tooth", INST),
             ("beside", "beside a tooth", AMB),
             ("control", "control position", SIG)]
    amps = sorted({k[1] for k in tally})

    fig, ax = plt.subplots(figsize=(8.6, 3.0))
    width = 0.26
    for i, (kind, label, col) in enumerate(kinds):
        vals, ns = [], []
        for a in amps:
            killed, total = tally[(kind, a)]
            vals.append(100.0 * killed / total if total else 0.0)
            ns.append(total)
        xs = [j + (i - 1) * width for j in range(len(amps))]
        ax.bar(xs, vals, width * 0.92, color=col, alpha=0.85,
               label=f"{label}  (n={sum(ns):,})")
        for x, v in zip(xs, vals):
            ax.text(x, v + 0.07, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=7, color=DIM, family="monospace")
    ax.set_xticks(range(len(amps)))
    ax.set_xticklabels([f"{a:g} mag" for a in amps])
    ax.set_xlabel("injected amplitude")
    ax.set_ylabel("false-kill rate  (%)")
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.5)
    ax.set_title("Real signals wrongly killed, 17,372 injections", loc="left",
                 fontsize=9.5, color=INK, pad=10)
    save(fig, "fig3_injection.png")


def fig_validity():
    """Where the method stops being trustworthy: period vs sector baseline."""
    rows = list(csv.DictReader(open(VAL / "v6_validity.csv")))
    mid = [(float(r["ratio_lo"]) + float(r["ratio_hi"])) / 2 for r in rows]
    rec = [float(r["median_recovery"]) for r in rows]
    fk = [float(r["false_kill_pct"]) for r in rows]

    fig, ax = plt.subplots(figsize=(8.6, 3.2))
    ax.axvspan(0.45, max(mid) + 0.05, color=INST, alpha=0.07, zorder=0)
    ax.plot(mid, rec, "-o", color=SIG, lw=2, ms=5, label="median signal recovery")
    ax.set_ylim(0.62, 1.02)
    ax.set_xlim(0, max(mid) + 0.05)
    ax.set_xlabel("period / sector baseline")
    ax.set_ylabel("recovery fraction", color=SIG)
    ax.tick_params(axis="y", colors=SIG)
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    ax2 = ax.twinx()
    ax2.plot(mid, fk, "--s", color=INST, lw=1.6, ms=4, label="false-kill rate")
    ax2.set_ylabel("false-kill (%)", color=INST)
    ax2.tick_params(axis="y", colors=INST)
    ax2.set_ylim(0, 9)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_color(INST)

    ax.text(0.62, 0.655, "outside validity domain", fontsize=8.5, color=INST,
            ha="center")
    # inline labels instead of a legend box, which collides with the curves
    ax.annotate("median signal recovery", (mid[1], rec[1]), xytext=(6, -16),
                textcoords="offset points", fontsize=8.5, color=SIG)
    ax2.annotate("false-kill rate", (mid[-2], fk[-2]), xytext=(-8, 10),
                 textcoords="offset points", fontsize=8.5, color=INST, ha="right")
    ax.set_title("Validity domain: recovery degrades as P approaches the sector length",
                 loc="left", fontsize=9.5, color=INK, pad=10)
    save(fig, "fig4_validity.png")


if __name__ == "__main__":
    print("writing figures:")
    fig_comb()
    fig_plane()
    fig_injection()
    fig_validity()

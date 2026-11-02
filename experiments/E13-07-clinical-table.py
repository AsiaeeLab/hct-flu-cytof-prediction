#!/usr/bin/env python3
"""Supplementary table of the clinical score statistics (October 2026 round).

Reads the E13-02 summaries and the E13-04 training-aware null p-values and
writes paper/tables/clinical_score.tex. Every number comes from those CSVs.
"""
import csv, os

CS = os.path.join(os.environ.get("FLU_RESULTS_ROOT",
                                 "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"),
                  "post_hoc_5000/clinical_score")
OUT = "/home/amir/repo/flu-cellular-prediction/paper/tables/clinical_score.tex"


SUBSETS = [("All", "All recipients"), ("Children", "Children"), ("Adults", "Adults")]


def read(name):
    with open(os.path.join(CS, name)) as f:
        return list(csv.DictReader(f))


def num(x, d=2):
    v = float(x)
    s = f"{abs(v):.{d}f}"
    return ("$-$" if v < 0 and float(s) != 0 else "") + s


def pv(x):
    v = float(x)
    if v < 0.001:
        return "$<$.001"
    if v <= 0.0010001:
        return ".001"
    out = f"{v:.3f}" if v < 0.01 else f"{v:.2f}" if v >= 0.1 else f"{v:.3f}"
    return out[1:] if out.startswith("0.") else out


null = {(r["statistic"], r["subset"]): r for r in read("null_pvalues.csv")}
def pn(stat, sub="All"):
    r = null.get((stat, sub))
    if r is None:
        return "--"
    p = float(r["p_null"])
    n = int(r["n_null"])
    return f"$\\le${1/(n+1):.3f}".replace("0.", ".", 1) if p <= 1 / (n + 1) + 1e-12 else pv(p)

rows = []
bt = {(r["subset"], r["quantity"]): r for r in read("bottom_fifth.csv")}
for sub, lab in [("All", "All recipients"), ("Children", "Children"), ("Adults", "Adults")]:
    lo, re_, df = bt[(sub, "rate_bottom_fifth")], bt[(sub, "rate_rest")], bt[(sub, "difference")]
    rows.append(("Fourfold rise, lowest fifth vs other fifths (\\%)", lab,
                 f"{100*float(lo['estimate']):.0f} vs {100*float(re_['estimate']):.0f}; difference {100*float(df['estimate']):.0f} [{100*float(df['ci_low']):.0f}, {100*float(df['ci_high']):.0f}]",
                 "--", pn("bottom_fifth_difference", sub) if sub == "All" else "--"))
lg = read("logistic_models.csv")
fe = read("logistic_first_episode.csv")[0]
# When every episode is a first episode, that subset repeats the all-episode row.
FIRST_EPISODE_SUBSET_DIFFERS = False
def lrow(r, lab, stat_sub):
    return ("Odds ratio for fourfold rise per SD", lab,
            f"{float(r['or']):.2f} [{float(r['or_lo']):.2f}, {float(r['or_hi']):.2f}]", pv(r["p"]),
            pn("logistic_or_per_sd", stat_sub) if stat_sub else "--")
m = {r["model"] + "|" + r["term"]: r for r in lg}
rows.append(lrow(m["baseline + score|z_comp"], "All recipients", "All"))
rows.append(lrow(m["baseline + score, Children|z_comp"], "Children", "Children"))
rows.append(lrow(m["baseline + score, Adults|z_comp"], "Adults", "Adults"))
# In the first-season analysis every episode is a first episode, so the
# "first episodes" subset repeats the "all episodes" row and is left out.
if FIRST_EPISODE_SUBSET_DIFFERS:
    rows.append(lrow(fe, "First episodes", "First episodes"))
rows.append(lrow(m["baseline + score + months + dose + cohort + age|z_comp"], "All, adjusted for clinical factors", None))
ia = m["baseline + score x cohort|z_comp:cohortAdults"]
rows.append(("Ratio of odds ratios, adults vs children", "Score-by-cohort interaction",
             f"{float(ia['or']):.2f} [{float(ia['or_lo']):.2f}, {float(ia['or_hi']):.2f}]", pv(ia["p"]), "--"))
fr = read("foldrise_per_sd.csv")[0]
rows.append(("Fold-rise ratio per SD at fixed baseline titer", "All recipients",
             f"{float(fr['ratio_per_sd']):.2f} [{float(fr['lo']):.2f}, {float(fr['hi']):.2f}]", pv(fr["p"]), "--"))
rc = {(r["subset"], r["statistic"]): r for r in read("rank_correlations.csv")}
for sub, lab in SUBSETS:
    r = rc[(sub, "Spearman")]
    rows.append(("Spearman correlation with composite fold-rise", lab,
                 f"{float(r['estimate']):.2f} [{float(r['ci_low']):.2f}, {float(r['ci_high']):.2f}]", "--", pn("spearman", sub)))
r = rc[("All", "Partial Spearman given baseline")]
rows.append(("Partial Spearman correlation given baseline titer", "All recipients",
             f"{float(r['estimate']):.2f} [{float(r['ci_low']):.2f}, {float(r['ci_high']):.2f}]", "--", pn("partial_spearman", "All")))
roc = {(r["subset"], r["quantity"]): r for r in read("roc.csv")}
for sub, lab in SUBSETS:
    a1, a3, d = roc[(sub, "AUC titers (M1)")], roc[(sub, "AUC titers + composition (M3)")], roc[(sub, "Delta AUC M1->M3")]
    rows.append(("AUC, titers (M1) vs titers + composition (M3)", lab,
                 f"{float(a1['estimate']):.2f} vs {float(a3['estimate']):.2f}; difference {float(d['estimate']):.3f} [{float(d['ci_low']):.3f}, {float(d['ci_high']):.3f}]",
                 pv(d["p_delong"]), pn("delta_auc_M1_M3", sub)))
a2, a4, d = roc[("All", "AUC titers + clinical (M2)")], roc[("All", "AUC titers + clinical + composition (M4)")], roc[("All", "Delta AUC M2->M4")]
rows.append(("AUC, titers + clinical (M2) vs all predictors (M4)", "All recipients",
             f"{float(a2['estimate']):.2f} vs {float(a4['estimate']):.2f}; difference {float(d['estimate']):.3f} [{float(d['ci_low']):.3f}, {float(d['ci_high']):.3f}]",
             pv(d["p_delong"]), "--"))
# composite outcome, pooled model evaluated within each cohort (E16-02)
with open(os.path.join(os.path.dirname(CS), "composite_main", "pooled_within_cohort.csv")) as f:
    cg = {(r["cohort"], r["comparison"]): dict(r, estimate=r["gain"], ci_low=r["lo"], ci_high=r["hi"], p_boot=r["p"])
          for r in csv.DictReader(f)}
for sub, lab in [("Children", "Children"), ("Adults", "Adults")]:
    for comp, clab in [("M1->M3", "M1$\\to$M3"), ("M2->M4", "M2$\\to$M4")]:
        r = cg[(sub, comp)]
        rows.append((f"Pooled model, MSE reduction {clab} (\\%)", lab,
                     f"{float(r['estimate']):.1f} [{float(r['ci_low']):.1f}, {float(r['ci_high']):.1f}]", pv(r["p_boot"]), "--"))

nnull = int(next(iter(null.values()))["n_null"]) if null else 0
body = "\n".join(f"{a} & {b} & {c} & {d} & {e} \\\\" for a, b, c, d, e in rows)
tex = r"""% Generated by experiments/E13-07-clinical-table.py
\begin{table}[p]
\centering
\caption{\textbf{Clinical summaries of the cellular composition score.}
Composite cellular scores are out-of-fold and standardized within
cross-validation fold; a fourfold rise is a geometric mean fold-rise across the
three antigens of at least four. Brackets give 95\% intervals: recipient-level
bootstrap ($B = 2{,}000$) for rates, correlations, AUCs and MSE reductions, and
recipient-clustered (sandwich) intervals for odds and fold-rise ratios.
Conventional $P$ values are Wald tests with recipient-clustered standard errors,
DeLong tests for AUCs, or recipient-bootstrap tests for MSE reductions.
Training-aware $P$ values come from """ + str(nnull) + r""" null-bootstrap replicates that repeat
resolution selection, penalty tuning and coefficient fitting under no cellular
effect beyond baseline titers; the smallest attainable value is
$1/(B+1)$. The odds ratios are adjusted for baseline titer; clinical factors are months since transplantation, vaccine dose, age and cohort.
MSE reductions are for the composite outcome.}
\label{tab:clinical_score}
\scriptsize
\renewcommand{\arraystretch}{0.95}
\begin{tabular}{p{4.7cm}p{2.5cm}p{4.5cm}p{1.0cm}p{1.3cm}}
\toprule
Statistic & Recipients & Estimate [95\% interval] & $P$ & Training-aware $P$ \\
\midrule
""" + body + r"""
\bottomrule
\end{tabular}
\end{table}
"""
open(OUT, "w").write(tex)
print("wrote", OUT, len(rows), "rows; null replicates", nnull)

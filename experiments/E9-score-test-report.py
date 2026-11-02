#!/usr/bin/env python3
"""Render all 30 exploratory score contrasts; requires all 999 null fits."""
from pathlib import Path
import math
import pandas as pd
root=Path('/home/amir/repo/flu-cellular-prediction')
import os
fast=Path(os.environ.get('FLU_RESULTS_ROOT','/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06'))
tests=pd.read_csv(fast/'score_null_bootstrap/score_tests.csv')
assert len(tests)==30 and set(tests.n_null)=={999}
contrasts=pd.read_csv(root/'paper/memo_figures/memo_tab1_score_contrasts.csv')
rows=[];report=[]
for virus,label in [('h1n1','H1N1'),('h3n2','H3N2'),('bvic','B/Victoria')]:
 for scope,stratum in [('Q1','Q1_low'),('Q2','Q2'),('Q3','Q3'),('Q4','Q4_high'),('Q1--Q3','Q1-Q3 pooled (within-quartile split)')]:
  a=tests[(tests.virus==virus)&(tests.scope==scope.replace('--','-'))].set_index('endpoint')
  c=contrasts[(contrasts.virus==virus)&(contrasts.stratum==stratum)].iloc[0]
  g=a.loc['log_fold_ratio'];r=a.loc['response_difference']
  assert abs(math.exp(g.estimate)-c.gmfr_ratio)<1e-10
  assert abs(r.estimate*100-c.seroconv_diff_pp)<1e-10
  row=f'{label} & {scope} & {c.gmfr_ratio:.2f} & {g.p_null:.3f} & {g.p_holm:.3f} & {c.seroconv_diff_pp:.1f} & {r.p_null:.3f} & {r.p_holm:.3f} '+r'\\'
  rows.append(row)
  report.append(dict(antigen=label,stratum=scope.replace('--','–'),fold_ratio=c.gmfr_ratio,fold_ratio_lo=c.gmfr_ratio_lo,fold_ratio_hi=c.gmfr_ratio_hi,p_fold=g.p_null,p_fold_holm=g.p_holm,rate_diff_pp=c.seroconv_diff_pp,rate_diff_lo=c.seroconv_diff_lo,rate_diff_hi=c.seroconv_diff_hi,p_rate=r.p_null,p_rate_holm=r.p_holm))
tex=r'''% Generated from timing-corrected null-bootstrap and contrast CSVs.
\begin{table}[p]
\centering
\caption{\textbf{Exploratory cellular-score contrasts and training-aware
null-bootstrap tests.} Fold ratio is the ratio of geometric mean observed
fold-rises (higher versus lower score). Rate difference is the corresponding
fourfold-response difference in percentage points. Each antigen uses its
own within-quartile score split; Q1--Q3 retains those assignments. Raw $p$
values compare observed separation with 999 recipient wild-bootstrap
replicates under a linear titers-only conditional-mean null, refitting
supervised resolution selection, M3 coefficients, penalty and score groups.
$p_{\rm H}$ is Holm-adjusted across all 30 comparisons. These working-null
tests incorporate training uncertainty; the conditional intervals in
Figures~\ref{fig:cellular_score} and~\ref{fig:score_supporting} do not.
Effects are descriptive and do not establish causal or individually
validated treatment-response differences.}
\label{tab:score_tests}
\footnotesize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{llrrrrrr}
\toprule
Antigen & Stratum & Fold ratio & $p$ & $p_{\rm H}$ & Rate diff. & $p$ & $p_{\rm H}$ \\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}
'''
(root/'paper/tables/score_tests.tex').write_text(tex)
pd.DataFrame(report).to_csv(root/'paper/review_2026-09-06/score_tests_readable.csv',index=False)
tests.to_csv(root/'paper/memo_figures/memo_tab7_score_tests.csv',index=False)
print(pd.DataFrame(report).to_string(index=False))

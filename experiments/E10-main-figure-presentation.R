#!/usr/bin/env Rscript
# Presentation-only revision: saved coordinates, estimates and intervals.
# No model fitting, resampling, atlas construction or LaTeX compilation.
# Run from repository root. Original figures remain available unchanged.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
fig_dir <- 'paper/images'
root <- FLU_ROOT
blue <- '#1B4F72'; orange <- '#E67E22'
th <- theme_minimal(base_size=8,base_family='sans') + theme(
  panel.grid.minor=element_blank(), axis.text=element_text(size=6.5,color='grey20'),
  axis.title=element_text(size=7.5), strip.text=element_text(face='bold',size=8),
  plot.title=element_text(face='bold',size=9), plot.tag=element_text(face='bold',size=11),
  legend.position='bottom', legend.title=element_blank())
savefig <- function(p,name,h) ggsave(file.path(fig_dir,paste0(name,'.png')),p,width=7.2,height=h,dpi=600,bg='white')
ann <- fread(FLU_ANN_FILE)
# Figure 1: schematic and cached UMAP, with original coordinates/assignments.
bio <- new.env(); load(FLU_BIO_FILE,bio)
u <- cbind(as.data.table(bio$dt_keep)[,.(meta=as.integer(meta))],as.data.table(bio$umap_dt))
K_POP <- max(u$meta); stopifnot(identical(sort(unique(u$meta)),seq_len(K_POP)))
FIRST_EPISODE <- grepl('first_episode|composite_', Sys.getenv('FLU_RESULTS_ROOT'))
COMPOSITE <- grepl('composite_', Sys.getenv('FLU_RESULTS_ROOT'))
u[,lineage:=ann$broad_lineage[match(meta,ann$cluster_id)]]
centers <- u[,.(UMAP1=median(UMAP1),UMAP2=median(UMAP2)),by=meta]
`%|%` <- function(a, b) ifelse(is.na(a), b, a)
lineage_colors <- setNames(c('#2166AC','#B2182B','#E6AB02','#1B9E77','#7570B3','#E7298A','#666666','#A6761D'),sort(unique(u$lineage)))
map <- ggplot(u,aes(UMAP1,UMAP2,color=lineage)) + geom_point(size=.08,alpha=.35) +
  geom_label_repel(data=centers,aes(UMAP1,UMAP2,label=meta),inherit.aes=FALSE,
    size=2.1,seed=19,label.padding=unit(.08,'lines'),min.segment.length=0,max.overlaps=Inf,
    box.padding=.25,point.padding=.05,segment.size=.2,fill='white') +
  scale_color_manual(values=lineage_colors,labels=function(x) c(`APC/Monocyte`='APC/monocyte',`NK cell`='NK-like',Other='Lineage\u2212',`Other T`='DNT')[x] %|% x) +
  guides(color=guide_legend(override.aes=list(size=2,alpha=1),nrow=1)) +
  labs(title='Immune-Cell Atlas',x='UMAP 1',y='UMAP 2') + th +
  theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
# Figure 1A: study design. Cohort composition, randomized dose arms, and the
# visit schedule. Numbers are from the episode-level trial metadata for the 343
# analyzed episodes. The V1-to-V2 interval is not drawn to scale because
# episode-level dosing dates were not available; V3 is four weeks after dose 2.
# V2 and V4 antibody titers exist in the trial data but are not used here, so
# they are drawn in grey: the panel states what was measured and what was modeled.
grey_line <- 'grey55'; grey_txt <- 'grey40'
cards <- data.table(
  xc    = c(26, 74),
  title = c('Pediatric cohort', 'Adult cohort'),
  l1    = c('148 recipients, 187 episodes', '115 recipients, 156 episodes'),
  l2    = c('Age 11.7 y (3.1–18.9)', 'Age 57.8 y (18.5–72.8)'),
  l3    = c('High dose 92 / standard dose 95', 'High dose 80 / standard dose 76'))
if (FIRST_EPISODE) cards <- data.table(
  xc    = c(26, 74),
  title = c('Pediatric cohort', 'Adult cohort'),
  l1    = c('147 recipients, first season', '113 recipients, first season'),
  l2    = c('Age 11.8 y (3.1–18.9)', 'Age 57.9 y (19.4–72.8)'),
  l3    = c('High dose 72 / standard dose 75', 'High dose 58 / standard dose 55'))
timing_label <- if (FIRST_EPISODE) 'Vaccination 3 to 36 months after transplantation (median 6); 74 second-season episodes held out for validation' else
  'Vaccination occurred 3 to 42 months after transplantation (median 9)'
visits <- data.table(
  x     = c(12, 37, 62, 88),
  tag   = c('V1', 'V2', 'V3', 'V4'),
  ev    = c('First dose', 'Second dose', 'Four weeks after dose 2', 'Later follow-up'),
  d1    = c('30-marker mass cytometry', 'HAI titers', 'HAI titers', 'HAI titers'),
  d2    = if (COMPOSITE) c('Baseline HAI titers', 'One-dose response', 'Model outcome', 'not used here') else c('Baseline HAI titers', 'not used here', 'Model outcome', 'not used here'),
  col   = c(blue, grey_txt, orange, grey_txt),
  col1  = c(blue, grey_txt, grey_txt, grey_txt),
  face  = c('bold', 'plain', 'bold', 'plain'))
scheme <- ggplot() +
  geom_rect(data=cards, aes(xmin=xc-23, xmax=xc+23, ymin=79, ymax=100),
    fill='#F2F5F7', color='grey65', linewidth=.3) +
  geom_text(data=cards, aes(xc, 96.5, label=title), size=2.7, fontface='bold') +
  geom_text(data=cards, aes(xc, 91.5, label=l1), size=2.35) +
  geom_text(data=cards, aes(xc, 87, label=l2), size=2.35) +
  geom_text(data=cards, aes(xc, 82.5, label=l3), size=2.35) +
  annotate('segment', x=c(26,74), xend=c(26,74), y=79, yend=76, linewidth=.35, color=grey_line) +
  annotate('segment', x=26, xend=74, y=76, yend=76, linewidth=.35, color=grey_line) +
  annotate('text', x=50, y=70,
    label=timing_label,
    size=2.35) +
  annotate('segment', x=12, xend=62, y=52, yend=52, linewidth=.4, color=grey_line) +
  annotate('segment', x=62, xend=88, y=52, yend=52, linewidth=.4, color=grey_line, linetype='22') +
  geom_point(data=visits, aes(x, 52, color=col), size=2.3) +
  geom_text(data=visits, aes(x, 57.5, label=tag, color=col), size=2.9, fontface='bold') +
  geom_text(data=visits, aes(x, 46.5, label=ev), size=2.3) +
  geom_text(data=visits, aes(x, 41.5, label=d1, color=col1), size=2.2) +
  geom_text(data=visits, aes(x, 37, label=d2, color=col, fontface=face), size=2.2) +
  scale_color_identity() +
  annotate('segment', x=12, xend=12, y=33.5, yend=27, linewidth=.45, color=blue) +
  annotate('segment', x=12, xend=62, y=27, yend=27, linewidth=.45, color=blue) +
  annotate('segment', x=62, xend=62, y=27, yend=33.5,
    arrow=arrow(length=unit(.07,'inches')), linewidth=.45, color=blue) +
  annotate('text', x=37, y=17,
    label=if (COMPOSITE) 'Predict the mean HAI log fold-rise from V1 to V3 across H1N1, H3N2 and B/Victoria' else 'Predict antibody log fold-rise (log(V3/V1)) for H1N1, H3N2 and B/Victoria',
    size=2.6, color=blue, fontface='bold') +
  coord_cartesian(xlim=c(0,100), ylim=c(10,101), clip='off') + theme_void() +
  labs(title='Study Design') +
  theme(plot.title=element_text(face='bold',size=9), plot.margin=margin(8,10,4,10),
    legend.position='none')
savefig((scheme / map)+plot_layout(heights=c(1.4,2.5))+
  plot_annotation(tag_levels='A', tag_suffix='.'),'E10_study_atlas',7.0)

# Figure 2: all six saved primary comparisons, grouped by information added.
ci <- fread(file.path(root,'post_hoc_5000/recipient_bootstrap/E3_recipient_boot.csv'))[scope=='overall']
key <- data.table(comparison=c('M1_to_M3','M2_to_M4','M1_to_M2','M3_to_M4','M1_to_M4','M2_to_M3'),
 label=c('Add immune composition to titers\nM1 → M3','Add immune composition to titers + clinical factors\nM2 → M4',
 'Add clinical factors to titers\nM1 → M2','Add clinical factors to titers + immune composition\nM3 → M4',
 'Add clinical factors and immune composition to titers\nM1 → M4','Compare titers + clinical factors with titers + immune composition\nM2 → M3'))
cc <- merge(key,ci,by='comparison',sort=FALSE); stopifnot(nrow(cc)==6)
cc[,label:=factor(label,levels=rev(key$label))]
cc[,primary:=comparison %in% c('M1_to_M3','M2_to_M4')]
ladder <- ggplot(cc,aes(delta_pct,label,color=primary)) +
 geom_vline(xintercept=0,linetype='dashed',color='grey55',linewidth=.35) +
 geom_linerange(aes(xmin=ci_low_pct,xmax=ci_high_pct),linewidth=.7) + geom_point(size=2.5) +
 scale_color_manual(values=c('TRUE'=blue,'FALSE'='#777777'),guide='none') +
 labs(x='Reduction in mean squared prediction error (%)',y=NULL) + th + theme(panel.grid.major.y=element_blank())
savefig(ladder,'E10_primary_prediction',3.8)

# Figure 3: use the exact feature summaries from the established plotting script.
# Evaluate its data preparation only, never its plotting/output section.
src <- readLines('experiments/E2-fig-feature-importance.R')
eval(parse(text=src[seq_len(grep('^# Plot ---',src)[1]-1)]))
agg[,virus:=factor(virus,levels=c('h1n1','h3n2','bvic'))]
levels(agg$display_label) <- gsub('CD38hi','CD38-high',levels(agg$display_label),fixed=TRUE)
feature_plot <- function(d,title) ggplot(d,aes(signed_imp,display_label,fill=signed_imp>0)) +
 geom_col(width=.75) + geom_vline(xintercept=0,linewidth=.3) +
 facet_wrap(~virus,nrow=1,labeller=labeller(virus=c(h1n1='H1N1',h3n2='H3N2',bvic='B/Victoria'))) +
 scale_fill_manual(values=c('TRUE'=orange,'FALSE'=blue),breaks=c('TRUE','FALSE'),labels=c('Higher predicted fold-rise','Lower predicted fold-rise')) +
 labs(x='Standardized signed coefficient',y=NULL,title=title) + th +
 theme(panel.grid.major.y=element_blank(),legend.text=element_text(size=7),axis.text.y=element_text(size=6))
savefig((feature_plot(agg[is_cluster==FALSE],'Baseline and Clinical Predictors') /
 feature_plot(agg[is_cluster==TRUE],'Immune-Cell Populations')) +
 plot_layout(heights=c(1.1,3.5),guides='collect')+plot_annotation(tag_levels='A', tag_suffix='.') & theme(legend.position='bottom'),'E10_predictors',6.7)

# Figure 4: unchanged estimates and CIs, group counts visible, pooled group offset.
groups <- fread('paper/memo_figures/memo_tab1_score_groups.csv')
ql <- c('Q1_low','Q2','Q3','Q4_high','Q1-Q3 pooled (within-quartile split)')
d <- groups[virus=='h3n2' & stratum %in% ql]
d[,x:=c(1,2,3,4,5.6)[match(stratum,ql)]]
d[,group:=factor(group,levels=c('High','Low'))]
ns <- dcast(d,stratum~group,value.var='n')[match(ql,stratum)]
xl <- paste0(c('Q1\nLowest','Q2\n','Q3\n','Q4\nHighest','Q1–Q3\nCombined'),'\nn=',ns$High,'/',ns$Low)
base <- ggplot(d,aes(x,color=group,shape=group,group=group)) +
 scale_color_manual(values=c(High=blue,Low=orange),labels=c('Higher cellular score','Lower cellular score')) +
 scale_shape_manual(values=c(High=16,Low=17),labels=c('Higher cellular score','Lower cellular score')) +
 scale_x_continuous(breaks=c(1,2,3,4,5.6),labels=xl) +
 geom_vline(xintercept=4.8,color='grey80',linewidth=.3) + th +
 theme(panel.grid.major.x=element_blank()) + labs(x='Baseline antibody quartile')
dodge <- position_dodge(width=.45)
p1 <- base + geom_hline(yintercept=4,linetype='dashed',color='grey55',linewidth=.35) +
 geom_linerange(aes(ymin=gmfr_lo,ymax=gmfr_hi),position=dodge,linewidth=.6) +
 geom_point(aes(y=gmfr),position=dodge,size=2.3) +
 scale_y_continuous(trans='log2',breaks=c(1,2,4,8,16,32,64,128)) + labs(y='Geometric mean H3N2 fold-rise')
p2 <- base + geom_linerange(aes(ymin=100*seroconv_lo,ymax=100*seroconv_hi),position=dodge,linewidth=.6) +
 geom_point(aes(y=100*seroconv),position=dodge,size=2.3) +
 scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25)) + labs(y='Recipients with a fourfold H3N2 rise (%)')
savefig((p1|p2)+plot_layout(guides='collect')+plot_annotation(tag_levels='A', tag_suffix='.') & theme(legend.position='bottom'),'E10_cellular_score',3.8)

# Figure 5: original cohort A/B; C becomes a separate supplemental figure.
src <- readLines('experiments/E5-figure.R')
eval(parse(text=src[seq_len(grep('^# Combine:',src)[1]-1)]))
# Primary pooled estimates now have their own main figure.
p_imp$data <- p_imp$data[cohort!='Pooled']
p_imp <- p_imp + labs(y='Reduction in prediction error (%)') +
  facet_wrap(~comparison,nrow=1,labeller=labeller(comparison=c('M1 → M3'='Add composition to titers\nM1 → M3','M2 → M4'='Add composition to titers + clinical factors\nM2 → M4'))) + th
levels(p_sel$data$label) <- gsub('CD38hi','CD38-high',levels(p_sel$data$label),fixed=TRUE)
p_sel <- p_sel + th + labs(x='Selection frequency')
savefig((p_imp/p_sel)+plot_layout(heights=c(1,2))+plot_annotation(tag_levels='A', tag_suffix='.'),'E10_cohort_comparison',6.4)
savefig(p_strat+th+labs(y='Reduction in prediction error (%)'),'E10_cohort_quartiles',3.4)
message('Presentation figures saved. No analyses refitted; no manuscript compiled.')

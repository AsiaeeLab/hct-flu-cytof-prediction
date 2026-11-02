#!/usr/bin/env Rscript
# Editorial changes to supplementary figure labels only.
# All numerical summaries and UMAP/model results are reused; no analyses refitted.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
s <- readLines('experiments/E8-review-figures.R')
eval(parse(text=s[1:(grep('^# Main:',s)-1)]))
p_labels <- function(virus_name,endpoint_name,scopes,label) label
# Reuse the saved episode scores and group summaries with their original intervals.
a <- grep('^# Supplement: continuous',s); b <- grep('^# Supplement: label agreement',s)
z <- s[a:(b-1)]
z <- z[!grepl('^save_plot',z)]
eval(parse(text=z))
continuous <- continuous + labs(x = 'H3N2 cellular score (out-of-fold)')
combined <- (continuous / all_antigens) + plot_layout(heights = c(1.2, 1)) + plot_annotation(tag_levels = 'A', tag_suffix = '.')
save_plot(combined,'E11_score_supporting.png',6.7)

# Keep the existing schematic proportions but remove claims about observed stages.
s <- readLines('experiments/E4-cartoon.R')
eval(parse(text=s[1:(grep('^p <- ggplot',s)-1)]))
stage_labels$label <- c('Earlier recovery','Intermediate','Later recovery')
p <- ggplot() +
 geom_rect(data=rect_dt,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=category),color='white',linewidth=.3) +
 geom_segment(data=arrows_dt,aes(x=x1,xend=x2,y=y,yend=y),arrow=arrow(length=unit(.18,'cm'),type='closed'),color='grey40',linewidth=.5) +
 scale_x_continuous(breaks=c(1,2,3),labels=stage_labels$label,limits=c(.4,3.6)) +
 scale_y_continuous(breaks=NULL,limits=c(0,116),expand=expansion(mult=c(.01,.03))) +
 scale_fill_manual(values=cat_colors,name=NULL,breaks=cat_levels,
 labels=function(x) gsub('Immature / DN','Immature / double-negative',x,fixed=TRUE)) +
 labs(x=NULL,y='Illustrative relative composition') + theme_minimal(base_size=9) +
 theme(panel.grid=element_blank(),legend.position='bottom',legend.text=element_text(size=7),axis.text.x=element_text(size=8)) +
 guides(fill=guide_legend(nrow=3,byrow=TRUE))
ggsave('paper/images/E11_reconstitution_schematic.png',p,width=7.2,height=4.3,dpi=600,bg='white')

REPEATER_WS <- file.path(FLU_ROOT, "post_hoc_5000/repeater_analysis/E4_repeater_workspace.RData")
if (file.exists(REPEATER_WS)) {
# Relabel the three-antigen PCA view using the saved coordinates.
e <- new.env(); load(file.path(FLU_ROOT, "post_hoc_5000/repeater_analysis/E4_repeater_workspace.RData"),e)
scores <- as.data.table(e$pca$x[,1:2]); scores[,obs_id:=rownames(e$pca$x)]
scores[,base_id:=sub('_yr[0-9]+$','',obs_id)]
scores[,year_abs:=as.integer(sub('.*_yr','',obs_id))]
scores[,returning:=base_id %in% e$rep_wide$base_id]
scores[,episode:=frank(year_abs,ties.method='first'),by=base_id]
scores <- merge(scores,e$rep_wide[,.(base_id,trans_h1n1,trans_h3n2,trans_bvic)],by='base_id',all.x=TRUE)
fg <- melt(scores[returning==TRUE],id.vars=c('base_id','PC1','PC2','episode'),
 measure.vars=c('trans_h1n1','trans_h3n2','trans_bvic'),variable.name='antigen',value.name='response')
fg[,antigen:=factor(antigen,levels=c('trans_h1n1','trans_h3n2','trans_bvic'),labels=c('H1N1','H3N2','B/Victoria'))]
trans_colors <- c('NR -> NR'='#d6604d','NR -> R'='#2ca02c','R -> NR'='#762a83','R -> R'='#1f77b4')
p <- ggplot() + geom_point(data=scores[returning==FALSE],aes(PC1,PC2),color='grey85',size=.5) +
 geom_path(data=fg[order(base_id,episode)],aes(PC1,PC2,group=base_id,color=response),
  arrow=arrow(length=unit(.10,'cm'),type='closed'),linewidth=.35,alpha=.8) +
 geom_point(data=fg,aes(PC1,PC2,color=response,shape=factor(episode)),size=1) +
 facet_wrap(~antigen,nrow=1) +
 scale_color_manual(values=trans_colors,name='Response trajectory',labels=function(x)gsub(' -> ',' → ',x)) +
 scale_shape_manual(values=c('1'=16,'2'=17),name='Episode') +
 labs(x=sprintf('PC1 (%.1f%%)',e$expl[1]),y=sprintf('PC2 (%.1f%%)',e$expl[2])) +
 theme_minimal(base_size=8) + theme(panel.grid.minor=element_blank(),legend.position='bottom',
 strip.text=element_text(face='bold'),legend.text=element_text(size=6.5),legend.title=element_text(size=7))
ggsave('paper/images/E11_repeater_pca.png',p,width=7.2,height=3.5,dpi=600,bg='white')

}
# Heatmap labels only: use cached scaled marker medians; no cell-level recomputation.
library(RColorBrewer)
med <- fread('paper/images/figS4_scaled_marker_medians.csv')
marker_names <- setdiff(names(med),'meta')
mat <- as.matrix(med[,..marker_names]);rownames(mat)<-paste0('C',med$meta)
annotation_dt <- fread(FLU_ANN_FILE)
annotation_dt[,annotation:=gsub('CD38hi','CD38-high',annotation,fixed=TRUE)]
s <- readLines('experiments/E7-heatmap-catalyst-scaling.R')
a<-grep('^row_hc <-',s);b<-grep('^ggsave',s)
eval(parse(text=s[a:(b-1)]))
ggsave('paper/images/E11_marker_heatmap.png',p,width=7.2,height=6.6,dpi=600,bg='white')

# Remove overlapping median annotation, retaining the reference line and estimates.
s <- readLines('experiments/E2-fig1-k-selection.R')
eval(parse(text=s[1:(grep('^# --- Combine panels',s)-1)]))
p_b$layers <- p_b$layers[-4]
# Panel A: the selection criterion is the mean over inner folds of the composite ratio (E16-01).
p_a$scales$scales <- Filter(function(x) !('fill' %in% x$aesthetics), p_a$scales$scales)
p_a <- p_a + scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b', midpoint = 1.0, name = 'Inner-fold\nMSE ratio')
# Panel C: the composite outcome (primary) replaces the mean over the three antigens.
ct <- as.data.table(readRDS(file.path(FLU_ROOT, 'post_hoc_5000/composite/composite_episode_table.rds')))
comp_perf <- as.data.frame(ct[, .(virus = 'Composite', mse_m1 = mean((y_comp - pred_M1)^2), mse_m3 = mean((y_comp - pred_M3)^2)),
                              by = .(fold = paste0('Fold ', fold))])
comp_perf$improvement <- 100 * (1 - comp_perf$mse_m3 / comp_perf$mse_m1)
perf_c <- rbind(as.data.frame(perf_df)[perf_df$virus != 'Overall', names(comp_perf)], comp_perf)
perf_c$virus <- factor(perf_c$virus, levels = c('H1N1', 'H3N2', 'B/Vic', 'Composite'))
p_c <- ggplot(perf_c, aes(x = fold, y = improvement, fill = virus)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) + geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c('H1N1' = '#66c2a5', 'H3N2' = '#fc8d62', 'B/Vic' = '#8da0cb', 'Composite' = '#1B4F72'),
                    labels = c('H1N1' = 'H1N1', 'H3N2' = 'H3N2', 'B/Vic' = 'B/Victoria', 'Composite' = 'Composite'), name = 'Outcome') +
  labs(x = NULL, y = expression('M1' %->% 'M3 error reduction (%)')) + theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = 'bottom')
print(perf_c[perf_c$virus == 'Composite', c('fold', 'improvement')])
fig1 <- p_a / (p_b|p_c) + plot_layout(heights=c(1,1.2)) + plot_annotation(tag_levels='A',tag_suffix='.') & theme(plot.tag = element_text(face = 'bold'))
ggsave('paper/images/E11_k_selection.png',fig1,width=9,height=6,dpi=600,bg='white')

if (file.exists(REPEATER_WS)) {
# Longitudinal labels use episode and the same CD38-high spelling as the paper.
s <- readLines('experiments/E8-review-figures.R')
a<-grep('^# Refresh the longitudinal',s);b<-grep('^ggsave\\(file.path\\(out, "E4_repeater_cluster',s)
eval(parse(text=s[a:(b-1)]))
delta_plot$data[,label:=gsub('CD38hi','CD38-high',label,fixed=TRUE)]
delta_plot <- delta_plot + labs(x='Mean CLR change (episode 2 − episode 1)')
ggsave('paper/images/E11_repeater_delta.png',delta_plot,width=9,height=8,dpi=300,bg='white')
counts <- repeater_env$rep_wide[, .N, by=.(n_responded_1,n_responded_2)]
p <- ggplot(counts,aes(factor(n_responded_1),factor(n_responded_2),fill=N)) +
 geom_tile() + geom_text(aes(label=N),color='white',size=4,fontface='bold') +
 scale_fill_gradient(low='#deebf7',high='#08306b') +
 labs(x='Episode 1: antigens with a fourfold rise',y='Episode 2: antigens with a fourfold rise') +
 theme_bw(base_size=10) + theme(legend.position='none')
ggsave('paper/images/E11_repeater_count.png',p,width=6,height=5,dpi=300,bg='white')
}

#!/usr/bin/env Rscript
# E14-09: Supplementary Table S2 (population annotations) for the first-season
# shared atlas. Labels are Josh Simmons's October 30, 2026 review, read from
# FLU_ANN_FILE; key markers follow the earlier curated lists (E14-03b). For each population carrying a label
# from a population with curated key markers, the same markers are listed with
# values recomputed on this atlas (mean expression, original scale); text
# qualifiers such as "CD27 low" are kept. Populations split from another
# population (", subset 2") list the three markers with the highest expression
# relative to the other populations (z-score across populations).
# Writes paper/tables/cluster_annotations.tex.
suppressPackageStartupMessages(library(data.table))
e <- new.env(); load(Sys.getenv("FLU_BIO_FILE"), e)
m <- as.data.frame(e$means_original); rownames(m) <- m$meta
M <- as.matrix(m[order(m$meta), e$marker_names]); Z <- scale(M)
ann <- fread(Sys.getenv("FLU_ANN_FILE"))
tex_old <- readLines("paper/archive/2026-10-14_pre-green-round/manuscript.tex")
i0 <- grep("label\\{tab:cluster_annotations\\}", tex_old); i1 <- i0 + grep("end\\{tabular\\}", tex_old[i0:length(tex_old)])[1]
rows <- grep("^\\s*[0-9]+\\s*&", tex_old[i0:i1], value = TRUE)
old <- rbindlist(lapply(strsplit(sub("\\\\\\\\\\s*$", "", rows), "&"), function(x)
  data.table(old_id = as.integer(trimws(x[1])), label = trimws(x[3]), mtxt = trimws(x[4]))))
fmt <- function(v) if (round(v, 2) >= 1) sprintf("%.1f", v) else sprintf("%.2f", v)
mk <- function(id, old_id, subset) {
  if (!subset) {
    key <- old$mtxt[match(old_id, old$old_id)]
    toks <- trimws(strsplit(key, ",")[[1]])
    return(paste(vapply(toks, function(t) if (grepl("=", t)) {
      nm <- sub("=.*", "", t); col <- nm
      if (col %in% colnames(M)) paste0(nm, "=", fmt(M[id, col])) else t } else t, ""), collapse = ", "))
  }
  top <- order(-Z[id, ])[1:3]
  paste(sprintf("%s=%s", colnames(M)[top], vapply(M[id, top], fmt, "")), collapse = ", ")
}
lab_tex <- function(s) { s <- gsub("\u2212", "$^-$", s, fixed = TRUE); gsub("+", "$^+$", s, fixed = TRUE) }
# Populations split from another population keep data-driven key markers.
pre <- fread(sub("cluster_annotations.csv", "cluster_annotations_pre_josh_2026-10-30.csv", Sys.getenv("FLU_ANN_FILE"), fixed = TRUE))
ann[, subset := grepl(", subset", pre$annotation[match(cluster_id, pre$cluster_id)])]
ann[, short := lab_tex(annotation)]
ann[, keys := mapply(mk, cluster_id, matched_old_cluster, subset)]
# Key markers from Josh Simmons's review (October 30, 2026), with values on this atlas.
# Low values mark markers the population lacks. Populations 3, 19 and 30 keep the
# curated text above.
josh_keys <- list(`1` = c("CD8", "CD45RO", "PD-1", "TIGIT"), `2` = c("CCR2", "HLA-DR", "CD69"),
  `4` = c("CD19", "CD21", "CD27"), `5` = c("CCR2", "HLA-DR", "CD69"), `6` = c("CD4", "CD45RO", "ICOS", "CD25"),
  `7` = c("CD4", "CD45RO", "CCR7", "CD27"), `8` = c("CD8", "CCR7", "CCR4", "CD45RO"), `9` = c("CD19", "OX40", "CD137", "CD138"),
  `10` = c("CD4", "CCR7", "CD25", "CD127"), `11` = c("HLA-DR", "CD86", "CD69"), `12` = c("CD19", "Ki-67"),
  `13` = c("CD19", "BTLA", "CXCR5", "CD21"), `14` = c("CD4", "CCR7", "CD45RO", "CD25"), `15` = c("CD8", "CD57", "CD27", "CD28"),
  `16` = c("CD4", "CD45RO", "CCR4", "CD25"), `17` = c("CCR2", "CD86", "Ki-67", "CD69"), `18` = c("CD8", "TIGIT", "CD45RO", "CCR7"),
  `20` = c("CD3", "CD4", "CD8", "CD38"), `21` = c("CD57", "CD38", "CD137", "CD3"), `22` = c("CD8", "TIGIT", "CD45RO", "CCR7"),
  `23` = c("HLA-DR", "CD86", "CCR2"), `24` = c("TIGIT", "CD38"), `25` = c("CD4", "CD45RO", "PD-1", "ICOS", "CXCR5"),
  `26` = c("CD19", "CXCR5", "CCR7"), `27` = c("CD8", "CCR7", "CD27", "CD45RO"), `28` = c("CD38", "Ki-67", "CD19", "CD20"),
  `29` = c("CD138", "CD38", "CXCR3"), `31` = c("CD19", "CD25", "CXCR5"), `32` = c("CD3", "CD38", "TIGIT", "Ki-67"),
  `33` = c("CD4", "CD57", "CD27"), `34` = c("CD19", "CCR7", "BTLA"), `35` = c("CCR2", "HLA-DR", "OX40", "CD137"))
stopifnot(all(unlist(josh_keys) %in% colnames(M)))
for (k in names(josh_keys)) { id <- as.integer(k); mk_j <- josh_keys[[k]]
  ann[cluster_id == id, keys := paste(sprintf("%s=%s", mk_j, vapply(M[id, mk_j], fmt, "")), collapse = ", ")] }
grp <- list(`B cell / Plasmablast` = c("B cell", "Plasmablast"), `CD4$^+$ T cell` = "CD4+ T", `CD8$^+$ T cell` = "CD8+ T",
            `APC / Monocyte` = "APC/Monocyte", `NK cell` = "NK cell", `Double-negative T / Lineage-negative` = c("Other T", "Other"))
lin_tex <- c(`B cell` = "B cell", Plasmablast = "Plasmablast", `CD4+ T` = "CD4$^+$ T", `CD8+ T` = "CD8$^+$ T",
             `APC/Monocyte` = "APC/Mono", `NK cell` = "NK", Other = "Lineage$^-$", `Other T` = "DN T")
body <- character()
for (g in names(grp)) {
  if (length(body)) body <- c(body, "    \\midrule")
  body <- c(body, sprintf("    \\multicolumn{4}{l}{\\textit{%s}} \\\\", g))
  x <- ann[broad_lineage %in% grp[[g]]][order(match(broad_lineage, grp[[g]]), cluster_id)]
  body <- c(body, sprintf("    %d & %s & %s & %s \\\\", x$cluster_id, lin_tex[x$broad_lineage], x$short, gsub("\\$\\^", "$^", x$keys)))
}
out <- c("\\begin{table}[p]", "  \\centering", "  \\color{revgreen}",
  "  \\captionsetup{font={small,color=revgreen},labelfont={bf,color=revgreen}}",
  "  \\caption{\\textbf{Cell-population annotations.}",
  sprintf("  Population number, broad lineage, cell type label, and key differentiating
  markers (mean expression on the original scale) for the %d populations of the
  shared atlas. Annotations were assigned using canonical marker profiles; the
  CyTOF panel lacks dedicated myeloid markers (CD14, CD16), so APC/monocyte
  annotations were inferred by exclusion (CD3$^-$CD19$^-$HLA-DR$^+$). Key markers
  distinguish a population from others of its lineage; low values mark markers it
  lacks. Low/high in labels refer to relative marker expression, not population
  frequency; double-negative
  T cells (DNT) have a CD3$^+$CD4$^-$CD8$^-$ phenotype. Plasmablasts are listed
  with B cells. CM, central memory; EM, effector memory; EMRA, effector
  memory re-expressing CD45RA, assigned here as CD45RO$^-$CCR7$^-$ because the
  panel has no CD45RA; Mono, monocyte; APC, antigen-presenting cell.}", nrow(ann)),
  "  \\label{tab:cluster_annotations}", "  \\scriptsize", "  \\begin{tabular}{@{}clll@{}}", "    \\toprule",
  "    \\# & Lineage & Cell type label & Key markers \\\\", "    \\midrule", body, "    \\bottomrule", "  \\end{tabular}", "\\end{table}")
writeLines(out, "paper/tables/cluster_annotations.tex")
cat(body, sep = "\n")

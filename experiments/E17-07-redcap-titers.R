#!/usr/bin/env Rscript
# E17-07: HAI titers from the REDCap exports (adult and pediatric HCT trials, HAI
# instrument only, as sent by Josh Simmons on October 29, 2026). This table
# replaces the titer columns of both older files: the May 2024 counts file codes
# a missing visit as 1 and had visit 4 copied into visit 1 for one second-season episode,
# and the July 2025 metadata blanked some B/Victoria results of <10.
#
# Each visit has two replicate titers per antigen. The episode titer is their
# geometric mean, with <10 recorded as 5, >=10240 as 10240, and NR or blank as
# missing. IVIG visit rows are not patient visits and are ignored. Study year:
# the adult trial began in 2017-18 and the pediatric trial in 2016-17 (year 1).
#
# Output: <raw>/hct_hai_redcap_2026-10-29.csv with id_yr, pid, yr and
# geo_mean_<antigen>_v<visit> for H1N1, H3N2, B/Victoria and B/Yamagata, visits 1-4.
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("methods", "00-paths.R"))
SRC <- file.path(paths$raw, "redcap_2026-10-29")
OUT <- file.path(paths$raw, "hct_hai_redcap_2026-10-29.csv")
SEASONS <- c("2016-2017", "2017-2018", "2018-2019", "2019-2020")
AG <- c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "BVIC", byam = "BYAM")

titer <- function(x) {
  x <- trimws(x)
  out <- suppressWarnings(as.numeric(x))
  out[x == "<10"] <- 5
  ge <- grepl("^>=", x); out[ge] <- as.numeric(sub("^>=", "", x[ge]))
  out[x %in% c("", "NR")] <- NA
  stopifnot(all(is.na(out) == (is.na(suppressWarnings(as.numeric(x))) & !(x == "<10" | ge))))
  out
}
read_redcap <- function(file, adult) {
  d <- fread(file.path(SRC, file), colClasses = "character", na.strings = NULL)
  d <- d[!grepl("IVIG", `which_visit?`)]
  stopifnot(all(grepl("^Visit [1-4]$", d$`which_visit?`)))
  num <- sub("^A?SCTF", "", d$`Study ID`)
  d[, pid := paste0(if (adult) "ASCT_" else "PSCT_", num)]
  d[, yr := match(`Year of Vaccine`, SEASONS) - if (adult) 1L else 0L]
  stopifnot(!anyNA(d$yr), all(d$yr >= 1))
  d[, visit := as.integer(sub("Visit ", "", `which_visit?`))]
  long <- rbindlist(lapply(names(AG), function(a) {
    r1 <- titer(d[[paste(AG[a], "result 1")]]); r2 <- titer(d[[paste(AG[a], "result 2")]])
    gm <- exp(rowMeans(log(cbind(r1, r2)), na.rm = TRUE)); gm[is.nan(gm)] <- NA
    data.table(pid = d$pid, yr = d$yr, visit = d$visit, antigen = a, gm = gm, n_rep = (!is.na(r1)) + (!is.na(r2)))
  }))
  stopifnot(!anyDuplicated(long[, .(pid, yr, visit, antigen)]))
  long
}
long <- rbind(read_redcap("adult_hct_hai_redcap.csv", TRUE), read_redcap("ped_hct_hai_redcap.csv", FALSE))
cat(sprintf("Visit-antigen titers: %d; single replicate: %d; missing: %d\n",
            nrow(long), sum(long$n_rep == 1), sum(is.na(long$gm))))
long[, col := sprintf("geo_mean_%s_v%d", antigen, visit)]
wide <- dcast(long, pid + yr ~ col, value.var = "gm")
cols <- sprintf("geo_mean_%s_v%d", rep(names(AG), 4), rep(1:4, each = 4))
for (c in setdiff(cols, names(wide))) wide[, (c) := NA_real_]
wide[, id_yr := paste0(pid, "_", yr)]
setcolorder(wide, c("id_yr", "pid", "yr", cols))
fwrite(wide[order(id_yr)], OUT)
cat(sprintf("Wrote %s: %d episodes\n", OUT, nrow(wide)))

# Agreement with the two older files on the analysis titers (visits 1-3).
cmp <- function(file, label) {
  old <- fread(file.path(paths$raw, file))
  use <- sprintf("geo_mean_%s_v%d", rep(c("h1n1", "h3n2", "bvic"), 3), rep(1:3, each = 3))
  m <- merge(wide[, c("id_yr", use), with = FALSE], old[, c("id_yr", use), with = FALSE], by = "id_yr", suffixes = c(".rc", ".old"))
  d <- rbindlist(lapply(use, function(c) {
    a <- m[[paste0(c, ".rc")]]; b <- m[[paste0(c, ".old")]]; b[!is.na(b) & b < 5] <- NA
    i <- which(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) > 1e-6))
    data.table(id_yr = m$id_yr[i], titer = rep(c, length(i)), redcap = a[i], old = b[i])
  }))
  cat(sprintf("%s: %d shared episodes, %d differing titers in %d episodes\n", label, nrow(m), nrow(d), uniqueN(d$id_yr)))
  print(d)
}
cmp("transplant_cytof_counts_hai_12MAY2024.csv", "May 2024 counts file (1 = missing)")
cmp("hsct_metadata_15JUL2025.csv", "July 2025 metadata")

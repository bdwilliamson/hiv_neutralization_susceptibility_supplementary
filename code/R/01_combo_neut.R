# set up datasets with individual and combination neutralization sensitivity
# for each bnAb combination of interest

# load required functions and packages -----------------------------------------
library("here")
library("readr")
library("janitor")
library("optparse")
library("glmnet")
library("data.table")
library("dplyr")
library("tidyr")

source(here::here("R", "00_utils.R"))

# define arguments -------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--bnab", default = "VRC07-523-LS+PGDM1400",
                     help = "The bnAb combination of interest")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

# load the datasets ------------------------------------------------------------
all_files <- list.files(here::here("..", "data"))
split_bnabs <- unlist(strsplit(args$bnab, "+", fixed = TRUE))
combo_bnab <- paste0(split_bnabs, collapse = "\\+")
has_combo_neut <- any(grepl(combo_bnab, all_files))
if (has_combo_neut) {
  bnabs <- c(split_bnabs, combo_bnab)
} else {
  bnabs <- split_bnabs
}
dat_list <- vector("list", length = length(bnabs))
for (i in seq_len(length(bnabs))) {
  bnab_file <- list.files(path = here::here("..", "data"),
                          pattern = paste0("_{1}", bnabs[i], "_{1}.*.csv"))
  dat <- read_csv(here::here("..", "data", bnab_file))
  dat$bnab <- janitor::make_clean_names(bnabs[i])
  dat_list[[i]] <- dat %>% 
    mutate(sens_ic50 = as.numeric(ic50 < 0),
           sens_ic80 = as.numeric(ic80 < 0)) %>% 
    select(-sens) %>% 
    group_by(seq.id.lanl) %>% 
    filter(!is.na(seq.id.catnap)) %>% 
    slice(1) %>% 
    ungroup()
  # some datasets may have duplicate rows?
}
long_dat <- fill_na_hxb2(data.table::rbindlist(dat_list, fill = TRUE))
long_x <- long_dat %>% 
  select(-ic50, -ic80, -iip, -starts_with("sens"), -bnab) 

# create a single dataset with individual and combined neutralization ----------
wide_dat <- long_dat %>% 
  pivot_wider(id_cols = seq.id.lanl, 
              values_from = matches("ic50|ic80|iip|sens_ic50|sens_ic80"),
              names_from = bnab) %>% 
  left_join(long_x %>%
              filter(complete.cases(long_x)) %>% 
              group_by(seq.id.lanl) %>% slice(1), 
            by = "seq.id.lanl")
nice_combo_bnab <- janitor::make_clean_names(paste0(split_bnabs, collapse = "_"))
nms_exclusions <- grepl("sens", names(wide_dat)) | grepl(nice_combo_bnab, names(wide_dat))
# note IC50, IC80 log10 transformed
all_ic50 <- 10 ^ as.matrix(wide_dat[, grepl("ic50_", names(wide_dat)) & !nms_exclusions]) 
all_ic80 <- 10 ^ as.matrix(wide_dat[, grepl("ic80_", names(wide_dat)) & !nms_exclusions])
wide_dat$ic50_additive <- log10(apply(all_ic50, 1, function(x) additive_combo_method(x)))
wide_dat$ic80_additive <- log10(apply(all_ic80, 1, function(x) additive_combo_method(x)))
wide_dat$iip_additive <- (-1) * log10(1 - get_iip(10 ^ wide_dat$ic50_additive, 10 ^ wide_dat$ic80_additive))
wide_dat$ic50_bh <- log10(predict_all_bh(conc = 0.5, ic50 = all_ic50, ic80 = all_ic80))
wide_dat$ic80_bh <- log10(predict_all_bh(conc = 0.8, ic50 = all_ic50, ic80 = all_ic80))
wide_dat$iip_bh <- (-1) * log10(1 - get_iip(10 ^ wide_dat$ic50_bh, 10 ^ wide_dat$ic80_bh))

wide_dat$sens_ic50_additive <- as.numeric(10 ^ wide_dat$ic50_additive < 1)
wide_dat$sens_ic80_additive <- as.numeric(10 ^ wide_dat$ic80_additive < 1)
wide_dat$sens_ic50_bh <- as.numeric(10 ^ wide_dat$ic50_bh < 1)
wide_dat$sens_ic80_bh <- as.numeric(10 ^ wide_dat$ic80_bh < 1)

wide_dat2 <- wide_dat %>% 
  select(starts_with("seq.id"), starts_with("ic"), starts_with("iip"), starts_with("sens"), everything())

saveRDS(wide_dat2, file = here::here("..", "data", paste0("analysis_dataset_", janitor::make_clean_names(args$bnab), ".rds")))
print("Data cleaning complete.")

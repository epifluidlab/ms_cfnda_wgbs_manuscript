library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(foreach)
library(doParallel)

setwd("/home/xfj2191/epifluidlab/projects/cfdna_ms_wgbs_pilot/pilot_75_samples/severity_plot_2")

name.order <- read.table("names_order.pdds_severe_vs_mild.txt", sep="\t", header=FALSE)
rownames(name.order) <- name.order[,2]

pdds.data <- read.table("survival_pdds_days.tsv", sep="\t", header=TRUE)
meta.data <- read.table("sample_meta_ms_wgbs_20191231.with_merged_names.v2.txt", sep="\t", header=FALSE)
rownames(meta.data) <- paste(meta.data[,1], meta.data[,9], sep=".")

d <- read.table("cpgs.pdds_severe_vs_mild.no_missing.merge300bp.add_value.methy.bed.gz",
                sep="\t", header=FALSE, na.strings="NaN", stringsAsFactors=FALSE)

pdds.data.inter <- read.table("pdds_over_time.tsv", sep="\t", header=TRUE)
rownames(pdds.data.inter) <- paste(pdds.data.inter[,1], pdds.data.inter[,2], sep=":")

### Build the main methylation matrix
methy.mat <- NULL
for(i in seq(7, ncol(d), 2)) {
  j <- i + 1
  s <- as.numeric(d[, i] / d[, j])
  methy.mat <- cbind(methy.mat, s)
}
methy.mat <- t(methy.mat)
methy.mat <- as.matrix(methy.mat)

colnames(methy.mat) <- d[,4]
rownames(methy.mat) <- name.order[,2]

# Filter columns with <10% missing
methy.mat <- methy.mat[, colSums(is.na(methy.mat)) < (nrow(methy.mat) * 0.1)]

sample_ids <- rownames(methy.mat)

### Build sample_meta
sample_meta <- NULL
sample_ids_clean <- NULL
for(sample_id in sample_ids){
  pid <- meta.data[sample_id, 2]
  # require at least 3 time points
  if(nrow(pdds.data[pdds.data[,1] %in% pid, ]) >= 3){
    tmp <- cbind(sample_id, pdds.data[pdds.data[,1] %in% pid, 1:3])
    m   <- meta.data[meta.data[,2] %in% pid, c(4,5,7)]
    tmp <- cbind(tmp, m)
    sample_meta <- rbind(sample_meta, tmp)
    sample_ids_clean <- c(sample_ids_clean, sample_id)
  }
}

methy.mat.clean <- methy.mat[sample_ids_clean, ]

## severity info for p-value calculation
patient_ids <- sample_meta[,2]
times       <- sample_meta[,4]
severity    <- sample_meta[,3]
gender      <- as.integer(as.factor(sample_meta[,5]))
race        <- as.integer(as.factor(sample_meta[,6]))
age         <- sample_meta[,7]

## For plotting: define time_points and so on
time_points <- seq(0,2000,50)
patient_ids_impute <- rep(unique(sample_ids_clean), each=length(time_points))
times_impute       <- rep(time_points, times=length(unique(sample_ids_clean)))

x <- cbind(patient_ids_impute, time_points,
           meta.data[patient_ids_impute,2],
           as.integer(as.factor(meta.data[patient_ids_impute,4])),
           as.integer(as.factor(meta.data[patient_ids_impute,5])),
           meta.data[patient_ids_impute,7])
x.names <- paste(x[,3], x[,2], sep=":")

severity_impute <- pdds.data.inter[x.names, 3]
gender_impute   <- as.integer(x[,4])
race_impute     <- as.integer(x[,5])
age_impute      <- as.integer(x[,6])

num_cols <- ncol(methy.mat.clean)

### Prepare parallel
#num_cores <- parallel::detectCores() - 1
num_cores=10
cl <- makeCluster(num_cores)
registerDoParallel(cl)

### We'll store results from each iteration
out_list <- foreach(i=1:num_cols, .combine='rbind',
                    .packages=c("lme4","lmerTest","dplyr","ggplot2")) %dopar% {
  # Print progress every 10,000 columns (may be interleaved in console)
  if(i %% 1000 == 0){
    cat("Processed column i =", i, "\n")
  }

  da <- methy.mat.clean[, i]
  threshold <- median(da, na.rm=TRUE)
  sample_ids_clean_local <- rownames(methy.mat.clean)  # to ensure scoping
  names(da) <- sample_ids_clean_local

  # build biomarker_values
  biomarker_values <- numeric(0)
  for(sample_id in sample_ids_clean_local){
    pid <- meta.data[sample_id, 2]
    n_rows_this <- nrow(sample_meta[sample_meta[,1]==sample_id,])
    biomarker_values <- c(biomarker_values, rep(da[sample_id], n_rows_this))
  }

  # define groups
  if(threshold == 0){
    groups <- ifelse(biomarker_values > threshold, 'High', 'Low')
  } else {
    groups <- ifelse(biomarker_values >= threshold, 'High', 'Low')
  }

  df_large <- data.frame(
    patient_id = patient_ids,
    time       = times,
    severity   = severity,
    biomarker_value = biomarker_values,
    gender     = gender,
    race       = race,
    age        = age,
    group      = groups
  )

  df_large <- df_large %>% filter(group != 'NA')
  if(length(unique(df_large$group)) == 1){
    # only one group => skip
    return(c(colname=colnames(methy.mat.clean)[i],
             base_e=NA, base_se=NA, base_p=NA, base_fdr=NA,
             interaction_e=NA, interaction_se=NA,
             interaction_p=NA, interaction_fdr=NA))
  }

  model <- lmer(severity ~ time * group + age + gender + race + (1 | patient_id), data=df_large)
  model_summary <- summary(model)
  interaction_term <- "time:groupLow"
  if(! (interaction_term %in% rownames(coef(model_summary)))) {
    return(c(colname=colnames(methy.mat.clean)[i],
             base_e=NA, base_se=NA, base_p=NA, base_fdr=NA,
             interaction_e=NA, interaction_se=NA,
             interaction_p=NA, interaction_fdr=NA))
  }

  interaction_p_value <- coef(model_summary)[interaction_term, "Pr(>|t|)"]
  fdr <- min(num_cols * interaction_p_value, 1)

  base_term <- "groupLow"
  base_e  <- if(base_term %in% rownames(coef(model_summary))) coef(model_summary)[base_term, "Estimate"] else NA
  base_se <- if(base_term %in% rownames(coef(model_summary))) coef(model_summary)[base_term, "Std. Error"] else NA
  base_p_value <- if(base_term %in% rownames(coef(model_summary))) coef(model_summary)[base_term, "Pr(>|t|)"] else NA
  base_fdr <- if(!is.na(base_p_value)) min(num_cols * base_p_value, 1) else NA

  if(fdr < 0.01){
    # build an impute data frame and save a PDF
    biomarker_values_impute <- rep(da, each=length(time_points))
    if(threshold==0){
      groups_impute <- ifelse(biomarker_values_impute > threshold, 'High', 'Low')
    } else {
      groups_impute <- ifelse(biomarker_values_impute >= threshold, 'High', 'Low')
    }

    df_large_impute <- data.frame(
      patient_id = patient_ids_impute,
      time       = times_impute,
      severity   = severity_impute,
      biomarker_value = biomarker_values_impute,
      gender     = gender_impute,
      race       = race_impute,
      age        = age_impute,
      group      = groups_impute
    ) %>% filter(group != 'NA')

    # create a pdf file
    file_name <- paste0("collection_prog_pdf.multithreads/sig_dmr_merge300bp.",
                        colnames(methy.mat.clean)[i],
                        ".threshold_", threshold, ".pdf")
    file_name <- gsub(":", "_", file_name)

    # Summaries
    plot_data_large <- df_large_impute %>%
      group_by(group, time) %>%
      summarise(
        mean_severity = mean(severity, na.rm=TRUE),
        sd_severity   = sd(severity, na.rm=TRUE),
        n = n(),
        se_severity   = sd_severity / sqrt(n)
      ) %>%
      ungroup() %>%
      mutate(
        ci_lower = mean_severity - 1.96 * se_severity,
        ci_upper = mean_severity + 1.96 * se_severity
      )

    x_position <- max(plot_data_large$time) * 0.05
    y_position <- max(plot_data_large$mean_severity, na.rm=TRUE) * 1.2

    # short function for star label
    format_p_value <- function(p_value){
      stars <- if(p_value < 0.001) "***" else if(p_value < 0.01) "**" else if(p_value < 0.05) "*" else "ns"
      paste0("P-value adj(Interaction)=", signif(p_value,3), " ", stars)
    }
    p_value_text <- format_p_value(fdr)

    pdf(file_name, paper="special", height=5, width=7)
    par(mar=c(5,5,1,6))

    p <- ggplot(plot_data_large, aes(x=time, y=mean_severity, color=group, fill=group)) +
      geom_line(size=1) +
      geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.2) +
      labs(x="Time", y="Mean PDDS", color="Methylation", fill="Methylation") +
      theme_minimal() +
      annotate("text", x=x_position, y=y_position, label=p_value_text, hjust=0, size=4)

    print(p)
    dev.off()
  }

  # return a row
  return(c(colname=colnames(methy.mat.clean)[i],
           base_e=base_e, base_se=base_se,
           base_p=base_p_value, base_fdr=base_fdr,
           interaction_e=coef(model_summary)[interaction_term, "Estimate"],
           interaction_se=coef(model_summary)[interaction_term, "Std. Error"],
           interaction_p=interaction_p_value,
           interaction_fdr=fdr))
}

# close cluster
stopCluster(cl)

# Convert out_list to data frame
out_df <- as.data.frame(out_list, stringsAsFactors=FALSE)
numeric_cols <- c("base_e","base_se","base_p","base_fdr","interaction_e","interaction_se","interaction_p","interaction_fdr")
for(nc in numeric_cols){
  out_df[[nc]] <- as.numeric(out_df[[nc]])
}

# Filter or store as needed
sig.ress <- subset(out_df, interaction_fdr < 0.01)

# Write results
write.table(sig.ress, "location_fdr001_model.correct_cov.multithreads.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Save environment
save.image("severity_plot.correct_cov.multithreads.RData")


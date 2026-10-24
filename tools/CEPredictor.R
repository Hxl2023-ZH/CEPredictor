#!/usr/bin/env Rscript

################################################################################
#                     CEPredictor : CE Potential Predictor
#                     Object-Oriented R6 + CLI Version
#
# Author : Xinliang Hu
# Version: 1.0.0
# Date   : 2026-10-22
#
# Usage:
#
# 1. Train model
# Rscript CEPredictor.R \
#   --input merged_annotation.txt \
#   --output train_results \
#   --train
#
# 2. Predict using pretrained model
# Rscript CEPredictor.R \
#   --input MAG_annotation_test.txt \
#   --output test_results \
#   --model train_result/models
#
################################################################################

suppressPackageStartupMessages({

  library(R6)
  library(optparse)

  library(data.table)
  library(tidyverse)

  library(randomForest)
  library(xgboost)
  library(caret)

  library(ggplot2)
  library(pROC)
  library(PRROC)

  library(FactoMineR)

})

set.seed(2026)

# rBOX is the mean of four module scores, three of which are 0/1 and one of
# which is 0/0.5/1, so it only takes multiples of 0.125. `>= 0.625` therefore
# is identical to `> 0.6` and includes MAGs sitting exactly on the boundary.
RBOX_CUTOFF <- 0.625

################################################################################
#                                Utility Functions
################################################################################

create_dir <- function(path){

  dir.create(
    path,
    recursive = TRUE,
    showWarnings = FALSE
  )

}

save_session_info <- function(outdir){

  writeLines(
    capture.output(sessionInfo()),
    file.path(outdir, "sessionInfo.txt")
  )

}

################################################################################
#                                Data Processor
################################################################################

DataProcessor <- R6Class(

  "DataProcessor",

  public = list(

    raw_df = NULL,
    gene_binary = NULL,
    MAGs = NULL,
    completeness = NULL,

    remove_cols = c(
      "H2ase",
      "BM",
      "ACOCT",
      "HadABC",
      "CarC",
      "HypCDEF"
    ),

    initialize = function(input_file){

      message("[INFO] Loading input data...")
      self$raw_df <- fread(input_file)

    },

    preprocess = function(){

      message("[INFO] Preprocessing data...")

      df <- self$raw_df %>%
        filter(Quality != "Low")

      remove_exist <- intersect(
        self$remove_cols,
        colnames(df)
      )

      if(length(remove_exist) > 0){

        df <- df %>%
          select(-all_of(remove_exist))

      }

      self$MAGs <- df$MAGs

      self$completeness <- df$completeness

      gene_table <- df %>%
        select(
          -MAGs,
          -completeness,
          -Quality
        )

      gene_binary <- gene_table
      gene_binary[gene_binary > 0] <- 1
      gene_binary[is.na(gene_binary)] <- 0
      gene_binary <- as.data.frame(gene_binary)

      rownames(gene_binary) <- self$MAGs

      self$gene_binary <- gene_binary

      message("[INFO] Data preprocessing completed.")

    }

  )

)

################################################################################
#                                Module Scorer
################################################################################

ModuleScorer <- R6Class(

  "ModuleScorer",

  public = list(

    module_list = NULL,
    module_matrix = NULL,

    initialize = function(){

      self$module_list <- list(

        condensation = c("ACAT"),
        reduction1 = c("HAD"),
        dehydration = c("ECH"),
        reduction2 = c("BCD","EtfAB"),
        product_release = c("ACT","CoAT"),
        lactate_pathway = c("LDH","LacP"),
        ethanol_pathway = c("ADH","ADA"),
        xylose_pathway = c("XI","XK","XylT")

      )

    },

    module_score = function(df, genes){

      valid <- intersect(
        genes,
        colnames(df)
      )

      if(length(valid) == 0){
        return(rep(0, nrow(df)))
      }

      mat <- as.matrix(
        df[, valid, drop = FALSE]
      )

      score <- rowSums(mat) / length(genes)

      return(score)

    },

    score = function(gene_binary){

      message("[INFO] Calculating module scores...")

      mat <- data.frame(
        row.names = rownames(gene_binary)
      )

      for(i in names(self$module_list)){

        mat[,i] <- self$module_score(
          gene_binary,
          self$module_list[[i]]
        )

      }

      mat$rBOX <- rowMeans(

        mat[, c(

          "condensation",
          "reduction1",
          "dehydration",
          "reduction2"

        )]

      )

      self$module_matrix <- mat

      message("[INFO] Module scoring completed.")

    }

  )

)

################################################################################
#                                CE Machine Learning Model
################################################################################

CEModel <- R6Class(

  "CEModel",

  public = list(

    xgb_model = NULL,
    rf_model = NULL,

    feature_names = NULL,

    core_genes = c(
      "ACAT",
      "HAD",
      "ECH",
      "BCD",
      "EtfAB"
    ),

    initialize = function(){
      message("[INFO] Initializing ML models...")
    },

    align_features = function(gene_binary){

      missing_cols <- setdiff(
        self$feature_names,
        colnames(gene_binary)
      )

      if(length(missing_cols) > 0){

        gene_binary[, missing_cols] <- 0

      }

      gene_binary <- gene_binary[, self$feature_names]

      return(gene_binary)

    },

    train = function(gene_binary, labels){

      message("[INFO] Training machine learning models...")

      feature_matrix <- gene_binary[, !colnames(gene_binary) %in%
                                      self$core_genes]

      self$feature_names <- colnames(feature_matrix)

      feature_matrix <- as.matrix(feature_matrix)

      set.seed(2026)

      train_index <- createDataPartition(

        labels,
        p = 0.8,
        list = FALSE

      )

      train_x <- feature_matrix[train_index,]
      train_y <- labels[train_index]

      test_x <- feature_matrix[-train_index,]
      test_y <- labels[-train_index]

      pos <- sum(train_y == 1)
      neg <- sum(train_y == 0)

      scale_pos_weight <- neg / pos

      message(
        paste(
          "[INFO] scale_pos_weight:",
          round(scale_pos_weight,3)
        )
      )

      # XGBoost
      dtrain <- xgb.DMatrix(

        train_x,
        label = train_y

      )

      dtest <- xgb.DMatrix(

        test_x,
        label = test_y

      )

      xgb_params <- list(

        objective = "binary:logistic",
        eval_metric = c("auc","aucpr"),
        max_depth = 5,
        eta = 0.1,
        subsample = 0.8,
        colsample_bytree = 0.8,
        scale_pos_weight = scale_pos_weight

      )

      xgb_args <- list(

        params = xgb_params,
        data = dtrain,
        nrounds = 500,
        early_stopping_rounds = 30,
        verbose = 0

      )

      # xgboost >= 2.0 renamed `watchlist` to `evals`
      if("evals" %in% names(formals(xgb.train))){

        xgb_args$evals <- list(
          train = dtrain,
          eval = dtest
        )

      }else{

        xgb_args$watchlist <- list(
          train = dtrain,
          eval = dtest
        )

      }

      self$xgb_model <- do.call(
        xgb.train,
        xgb_args
      )

      # RandomForest

      self$rf_model <- randomForest(

        x = train_x,
        y = as.factor(train_y),
        ntree = 1000,
        importance = TRUE,
        classwt = c(
          "0" = 1,
          "1" = scale_pos_weight
        )

      )

      message("[INFO] Model training completed.")

    },

    predict = function(gene_binary){

      message("[INFO] Predicting CE probabilities...")

      gene_binary <- self$align_features(
        gene_binary
      )

      feature_matrix <- gene_binary[, self$feature_names]
      feature_matrix <- as.matrix(feature_matrix)

      prob_xgb <- predict(
        self$xgb_model,
        feature_matrix
      )

      prob_rf <- predict(

        self$rf_model,
        feature_matrix,
        type = "prob"

      )[,2]

      ml_prob <- (
        prob_xgb + prob_rf
      ) / 2

      return(ml_prob)

    },

    save_model = function(model_dir = "models"){

      message("[INFO] Saving models...")

      create_dir(model_dir)

      # save xgb
      xgb.save(

        self$xgb_model,
        file.path(
          model_dir,
          "xgb_model.model"
        )

      )

      # save rf
      saveRDS(

        self$rf_model,
        file.path(
          model_dir,
          "rf_model.rds"
        )

      )

      # save metadata
      meta <- list(

        feature_names = self$feature_names,
        core_genes = self$core_genes,
        version = "CEPredictor v1.0.0",
        date = Sys.Date()

      )

      saveRDS(

        meta,

        file.path(
          model_dir,
          "model_meta.rds"
        )

      )

      message("[INFO] Models saved successfully.")

    },

    load_model = function(model_dir = "models"){

      message("[INFO] Loading pretrained models...")

      self$xgb_model <- xgb.load(

        file.path(
          model_dir,
          "xgb_model.model"
        )

      )

      self$rf_model <- readRDS(

        file.path(
          model_dir,
          "rf_model.rds"
        )

      )

      meta <- readRDS(

        file.path(
          model_dir,
          "model_meta.rds"
        )

      )

      self$feature_names <- meta$feature_names

      self$core_genes <- meta$core_genes

      message("[INFO] Models loaded successfully.")

    }

  )

)

################################################################################
#                                CEPS Calculator
################################################################################

CEPSCalculator <- R6Class(

  "CEPSCalculator",

  public = list(

    calculate = function(

      module_matrix,
      ml_prob,
      completeness

    ){

      message("[INFO] Calculating CEPS scores...")

      CEPS <- (

        0.5 * ml_prob +
        0.3 * module_matrix$product_release +
        0.2 * 0.01 * completeness

      )

      class <- rep("Low", length(CEPS))
      class[CEPS > 0.60] <- "Medium"
      class[CEPS > 0.80] <- "High"

      core_rule <- module_matrix$rBOX >= RBOX_CUTOFF

      result <- data.frame(

        MAG = rownames(module_matrix),
        CEPS = CEPS,
        class = class,
        core_rule = core_rule,
        ML_prob = ml_prob,
        completeness = completeness

      )

      return(result)

    }

  )

)

################################################################################
#                                Visualization
################################################################################

CEVisualizer <- R6Class(

  "CEVisualizer",

  public = list(

    outdir = NULL,
    initialize = function(outdir){
      self$outdir <- outdir
      create_dir(outdir)
    },

    plot_ceps_distribution = function(result_df){

      message("[INFO] Plotting CEPS distribution...")

      p <- ggplot(
        result_df,
        aes(
          x = class,
          y = CEPS,
          color = class
        )
      ) +

        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        theme_bw(base_size = 14) +

        scale_color_manual(
          values = c(
            Low = "gray",
            Medium = "#00688B",
            High = "#E64B35"
          )
        )

      ggsave(file.path(self$outdir, "Figure01_CEPS_distribution.pdf"), p, width = 5, height = 5)

    },

    plot_pca = function(gene_binary, result_df){

      message("[INFO] Plotting PCA...")

      pca <- PCA(gene_binary, graph = FALSE)
      coord <- as.data.frame(pca$ind$coord)
      coord$MAG <- rownames(coord)

      df <- merge(
        coord,
        result_df,
        by = "MAG"
      )

      p <- ggplot(df, aes(Dim.1, Dim.2, color = class)) +

        geom_point(alpha = 0.7) +
        theme_bw(base_size = 14) +

        scale_color_manual(
          values = c(
            Low = "gray",
            Medium = "#00688B",
            High = "#E64B35"
          )
        ) +

        xlab("PC1") +
        ylab("PC2")

      ggsave(file.path(self$outdir, "Figure02_PCA.pdf"), p, width = 6, height = 5)

    }

  )

)

################################################################################
#                                Main Pipeline
################################################################################

CEPipeline <- R6Class(

  "CEPipeline",

  public = list(

    input_file = NULL,
    output_dir = NULL,

    dp = NULL,
    scorer = NULL,
    model = NULL,
    calculator = NULL,
    viz = NULL,

    initialize = function(input_file, output_dir = "results"){

      self$input_file <- input_file
      self$output_dir <- output_dir

      create_dir(output_dir)
      create_dir(file.path(output_dir, "tables"))
      create_dir(file.path(output_dir, "figures"))

      self$dp <- DataProcessor$new(input_file)
      self$scorer <- ModuleScorer$new()
      self$model <- CEModel$new()
      self$calculator <- CEPSCalculator$new()
      self$viz <- CEVisualizer$new(

        file.path(
          output_dir,
          "figures"
        )

      )

    },

    ##########################################################################
    # TRAIN MODE

    run_train = function(){

      message("[INFO] Running training pipeline...")

      # preprocess
      self$dp$preprocess()

      # module score
      self$scorer$score(
        self$dp$gene_binary
      )

      # label generation (same rule as MAGs_CE_Potential_Classifer.R)
      labels <- ifelse(
        self$scorer$module_matrix$rBOX >= RBOX_CUTOFF,
        1,
        0
      )

      # model training
      self$model$train(
        self$dp$gene_binary,
        labels
      )

      # save model
      self$model$save_model(
        file.path(
          self$output_dir,
          "models"
        )
      )

      # prediction
      ml_prob <- self$model$predict(
        self$dp$gene_binary
      )

      # CEPS calculation
      result <- self$calculator$calculate(

        self$scorer$module_matrix,
        ml_prob,
        self$dp$completeness

      )

      # save results
      write.table(
        result,
        file.path(
          self$output_dir,
          "tables",
          "CEPS_prediction.txt"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      # visualization
      self$viz$plot_ceps_distribution(result)

      self$viz$plot_pca(
        self$dp$gene_binary,
        result
      )

      # session info
      save_session_info(
        self$output_dir
      )
      message("[INFO] Training pipeline completed.")
    },

    ##########################################################################
    # PREDICTION MODE

    run_predict = function(model_dir){

      message("[INFO] Running prediction pipeline...")

      # preprocess

      self$dp$preprocess()

      # module score
      self$scorer$score(self$dp$gene_binary)

      # load pretrained model
      self$model$load_model(model_dir)

      # predict
      ml_prob <- self$model$predict(self$dp$gene_binary)

      # calculate CEPS
      result <- self$calculator$calculate(self$scorer$module_matrix, ml_prob, self$dp$completeness)

      # save
      write.table(
        result,
        file.path(
          self$output_dir,
          "tables",
          "CEPS_prediction.txt"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      # plot
      self$viz$plot_ceps_distribution(
        result
      )
      self$viz$plot_pca(
        self$dp$gene_binary,
        result
      )

      # session info
      save_session_info(
        self$output_dir
      )
      message("[INFO] Prediction pipeline completed.")

    }

  )

)

################################################################################
#                                CLI Argument Parser
################################################################################

option_list <- list(

  make_option(

    c("-i","--input"),
    type = "character",
    help = "Input annotation table"

  ),

  make_option(

    c("-o","--output"),
    type = "character",
    default = "results",
    help = "Output directory"

  ),

  make_option(

    c("-m","--model"),
    type = "character",
    default = NULL,
    help = "Pretrained model directory"

  ),

  make_option(

    c("--train"),
    action = "store_true",
    default = FALSE,
    help = "Train new model"

  )

)

opt <- parse_args(OptionParser(option_list = option_list))

################################################################################
#                                Run Pipeline
################################################################################

pipeline <- CEPipeline$new(

  input_file = opt$input,
  output_dir = opt$output

)

###############################################################################
# Train mode
###############################################################################

if(opt$train){

  pipeline$run_train()

}else{

  if(is.null(opt$model)){
    stop("Please provide pretrained model directory using --model")
  }
  pipeline$run_predict(model_dir = opt$model)

}

message("[INFO] All analysis completed successfully.")
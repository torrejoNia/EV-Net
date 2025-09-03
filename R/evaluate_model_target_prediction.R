
#' @title Evaluation of target gene prediction.
#'
#' @description \code{evaluate_target_prediction} Evaluate how well the model (i.e. the inferred EV_cargo-target probability scores) is able to predict the observed response to a EV_cargo (e.g. the set of DE genes after treatment of cells by a EV_cargo). It shows several classification evaluation metrics for the prediction. Different classification metrics are calculated depending on whether the input EV_cargo-target matrix contains probability scores for targets or discrete target assignments.
#'
#' @usage
#' evaluate_target_prediction(setting,EV_cargo_target_matrix, EV_cargo_position = "cols")
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the EV_cargo(s) active in the setting of interest; .$response: named logical vector indicating whether a target is a TRUE target of the possibly active EV_cargo(s) or a FALSE.
#' @param EV_cargo_target_matrix A matrix of EV_cargo-target probabilty scores (or discrete target assignments).
#' @param EV_cargo_position Indicate whether the EV_cargo in the EV_cargo-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"

#' @return A data.frame with following variables: setting, EV_cargo nd for probabilistic predictions: auroc, aupr, aupr_corrected (aupr - aupr for random prediction), sensitivity_roc (proxy measure, inferred from ROC), specificity_roc (proxy measure, inferred from ROC), mean_rank_GST_log_pval (-log10 of p-value of mean-rank gene set test), pearson (correlation coefficient), spearman (correlation coefficient); whereas for categorical predictions: accuracy, recall, specificity, precision, F1, F0.5, F2, mcc, informedness, markedness, fisher_pval_log (which is -log10 of p-value fisher exact test), fisher odds.\cr
#' "mean_rank_GST_log_pval" will only be included in the dataframe if limma is installed. From NicheNet v2.1.7 onwards, limma is no longer a hard dependency of NicheNet.
#'
#'
#' @importFrom ROCR prediction performance
#' @importFrom caTools trapz
#' @importFrom data.table data.table
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' setting = lapply(expression_settings_validation[1],convert_expression_settings_evaluation)
#' EV_cargo = extract_EV_cargo_from_settings(setting)
#' EV_cargo_target_matrix = construct_EV_cargo_target_matrix(weighted_networks, EV_cargo)
#' perf1 = lapply(setting,evaluate_target_prediction,EV_cargo_target_matrix)
#' print(head(perf1))
#' perf2 = lapply(setting,evaluate_target_prediction,make_discrete_EV_cargo_target_matrix(EV_cargo_target_matrix))
#' }
#' @export
#'
evaluate_target_prediction = function(setting,EV_cargo_target_matrix, EV_cargo_position = "cols"){
  ## still make evaluation multiple EV_cargo possible
  # input check
  if (!is.list(setting))
    stop("setting must be a list")
  if(!is.character(setting$from) | !is.character(setting$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(setting$response) | is.null(names(setting$response)))
    stop("setting$response should be named logical vector containing class labels of the response that needs to be predicted ")
  if(!is.matrix(EV_cargo_target_matrix))
    stop("EV_cargo_target_matrix should be a matrix")
  if(!is.double(EV_cargo_target_matrix) & !is.logical(EV_cargo_target_matrix))
    stop("EV_cargo_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (EV_cargo_position != "cols" & EV_cargo_position != "rows")
    stop("EV_cargo_position must be 'cols' or 'rows'")

  requireNamespace("dplyr")

  if (length(setting$from) == 1){
    EV_cargo_oi = setting$from
  } else {
    EV_cargo_oi = paste0(setting$from,collapse = "-")
  }
  if (EV_cargo_position == "cols"){
    if((EV_cargo_oi %in% colnames(EV_cargo_target_matrix)) == FALSE)
      stop("EV_cargo should be in EV_cargo_target_matrix")
    prediction_vector = EV_cargo_target_matrix[,EV_cargo_oi]
    names(prediction_vector) = rownames(EV_cargo_target_matrix)
  } else if (EV_cargo_position == "rows") {
    if((EV_cargo_oi %in% rownames(EV_cargo_target_matrix)) == FALSE)
      stop("EV_cargo should be in EV_cargo_target_matrix")
    prediction_vector = EV_cargo_target_matrix[EV_cargo_oi,]
    names(prediction_vector) = colnames(EV_cargo_target_matrix)
  }
  response_vector = setting$response

  if(sd(prediction_vector) == 0)
    warning("all target gene probability score predictions have same value")
  if(sd(response_vector) == 0)
    stop("all genes have same response")
  performance = evaluate_target_prediction_strict(response_vector,prediction_vector,is.double(prediction_vector))
  output = bind_cols(tibble(setting = setting$name, EV_cargo = EV_cargo_oi), performance)

  return(output)
}

#' @title Evaluation of target gene prediction for multiple EV_cargo.
#'
#' @description \code{evaluate_multi_EV_cargo_target_prediction} Evaluate how well a trained model is able to predict the observed response to a combination of EV_cargo (e.g. the set of DE genes after treatment of cells by multiple EV_cargo). A classificiation algorithm chosen by the user is trained to construct one model based on the target gene predictions of all EV_cargo of interest (EV_cargo are considered as features). Several classification evaluation metrics for the prediction are calculated depending on whether the input EV_cargo-target matrix contains probability scores for targets or discrete target assignments. In addition, variable importance scores can be extracted to rank the possible active EV_cargo in order of importance for response prediction.
#'
#' @usage
#' evaluate_multi_EV_cargo_target_prediction(setting,EV_cargo_target_matrix, EV_cargo_position = "cols", algorithm, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4,ignore_errors = FALSE,continuous = TRUE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the EV_cargo(s) active in the setting of interest; .$response: named logical vector indicating whether a target is a TRUE target of the possibly active EV_cargo(s) or a FALSE.
#' @param EV_cargo_target_matrix A matrix of EV_cargo-target probabilty scores (recommended) or discrete target assignments (not-recommended).
#' @param EV_cargo_position Indicate whether the EV_cargo in the EV_cargo-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#' @param algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost", "pls" (load "pls" package before!); without: "lda","naive_bayes", "pcaNNet". Please notice that not all these algorithms work when the features (i.e. EV_cargo vectors) are categorical (i.e. discrete class assignments).
#' @param var_imps Indicate whether in addition to classification evaluation performances, variable importances should be calculated. Default: TRUE.
#' @param cv Indicate whether model training and hyperparameter optimization should be done via cross-validation. Default: TRUE. FALSE might be useful for applications only requiring variable importance, or when final model is not expected to be extremely overfit.
#' @param cv_number The number of folds for the cross-validation scheme: Default: 4; only relevant when cv == TRUE.
#' @param cv_repeats The number of repeats during cross-validation. Default: 2; only relevant when cv == TRUE.
#' @param parallel Indiciate whether the model training will occur parallelized. Default: FALSE. TRUE only possible for non-windows OS.
#' @param n_cores The number of cores used for parallelized model training via cross-validation. Default: 4. Only relevant on non-windows OS.
#' @param ignore_errors Indiciate whether errors during model training by caret should be ignored such that another model training try will be initiated until model is trained without raising errors. Default: FALSE.
#' @param continuous Indicate whether during training of the model, model training and evaluation should be done on class probabilities or discrete class labels. For huge class imbalance, we recommend setting this value to TRUE. Default: TRUE.
#'
#' @return A list with the following elements. $performances: data frame containing classification evaluation measure for classification on the test folds during training via cross-validation; $performances_training: data frame containing classification evaluation measures for classification of the final model (discrete class assignments) on the complete data set (performance can be severly optimistic due to overfitting!); $performance_training_continuous: data frame containing classification evaluation measures for classification of the final model (class probability scores) on the complete data set (performance can be severly optimistic due to overfitting!) $var_imps: data frame containing the variable importances of the different EV_cargo (embbed importance score for some classification algorithms, otherwise just the auroc); $prediction_response_df: data frame containing for each gene the EV_cargo-target predictions of the individual EV_cargo, the complete model and the response as well; $setting: name of the specific setting that needed to be evaluated; $EV_cargo: EV_cargo of interest.
#'
#' @importFrom ROCR prediction performance
#' @importFrom caTools trapz
#' @import caret
#' @importFrom purrr safely
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' setting = convert_expression_settings_evaluation(expression_settings_validation$TGFB_IL6_timeseries) %>% list()
#' EV_cargo = extract_EV_cargo_from_settings(setting)
#' EV_cargo_target_matrix = construct_EV_cargo_target_matrix(weighted_networks, EV_cargo)
#' output = lapply(setting,evaluate_multi_EV_cargo_target_prediction,EV_cargo_target_matrix,EV_cargo_position = "cols",algorithm = "glm")
#' output = lapply(setting,evaluate_multi_EV_cargo_target_prediction,make_discrete_EV_cargo_target_matrix(EV_cargo_target_matrix),EV_cargo_position = "cols",algorithm = "glm" )
#' }
#' @export
#'
evaluate_multi_EV_cargo_target_prediction = function(setting,EV_cargo_target_matrix, EV_cargo_position = "cols", algorithm, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE){
  if (!is.list(setting))
    stop("setting must be a list")
  if(!is.character(setting$from) | !is.character(setting$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(setting$response) | is.null(names(setting$response)))
    stop("setting$response should be named logical vector containing class labels of the response that needs to be predicted ")
  if(!is.matrix(EV_cargo_target_matrix))
    stop("EV_cargo_target_matrix should be a matrix")
  if(!is.double(EV_cargo_target_matrix) & !is.logical(EV_cargo_target_matrix))
    stop("EV_cargo_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (EV_cargo_position != "cols" & EV_cargo_position != "rows")
    stop("EV_cargo_position must be 'cols' or 'rows'")
  if(!is.character(algorithm))
    stop("algorithm should be a character vector")
  if(!is.logical(var_imps) | length(var_imps) > 1)
    stop("var_imps should be a logical vector: TRUE or FALSE")
  if(!is.logical(cv) | length(cv) > 1)
    stop("cv should be a logical vector: TRUE or FALSE")
  if(!is.numeric(cv_number) | length(cv_number) > 1)
    stop("cv_number should be a numeric vector of length 1")
  if(!is.numeric(cv_repeats) | length(cv_repeats) > 1)
    stop("cv_repeats should be a numeric vector of length 1")
  if(!is.logical(parallel) | length(parallel) > 1)
    stop("parallel should be a logical vector: TRUE or FALSE")
  if(!is.numeric(n_cores) | length(n_cores) > 1)
    stop("n_cores should be a numeric vector of length 1")
  if(!is.logical(ignore_errors) | length(ignore_errors) > 1)
    stop("ignore_errors should be a logical vector: TRUE or FALSE")
  if(!is.logical(continuous) | length(continuous) > 1)
    stop("continuous should be a logical vector: TRUE or FALSE")
  requireNamespace("dplyr")

  EV_cargo_oi = setting$from


  if (EV_cargo_position == "cols"){
    if(sum((EV_cargo_oi %in% colnames(EV_cargo_target_matrix)) == FALSE) > 0)
      stop("EV_cargo should be in EV_cargo_target_matrix")
    prediction_matrix = EV_cargo_target_matrix[,EV_cargo_oi]
    target_genes = rownames(EV_cargo_target_matrix)
  } else if (EV_cargo_position == "rows") {
    if(sum((EV_cargo_oi %in% rownames(EV_cargo_target_matrix)) == FALSE) > 0)
      stop("EV_cargo should be in EV_cargo_target_matrix")
    prediction_matrix = EV_cargo_target_matrix[EV_cargo_oi,] %>% t()
    target_genes = colnames(EV_cargo_target_matrix)
  }

  response_vector = setting$response
  if(sd(response_vector) == 0)
    stop("all genes have same response")
  response_df = tibble(gene = names(response_vector), response = response_vector %>% make.names() %>% as.factor())

  prediction_df = prediction_matrix %>% data.frame() %>% as_tibble()

  if(is.double(prediction_matrix) == FALSE){
    convert_categorical_factor = function(x){
      x = x %>% make.names() %>% as.factor()
    }
    prediction_df = prediction_df %>% mutate_all(funs(convert_categorical_factor))
  }

  prediction_df = tibble(gene = target_genes) %>% bind_cols(prediction_df)
  combined = inner_join(response_df,prediction_df, by = "gene")
  if (nrow(combined) == 0)
    stop("Gene names in response don't accord to gene names in EV_cargo-target matrix (did you consider differences human-mouse namings?)")

  train_data = combined %>% select(-gene) %>% rename(obs = response) %>% data.frame()

  output = wrapper_caret_classification(train_data,algorithm,continuous = continuous,var_imps,cv,cv_number,cv_repeats,parallel,n_cores,ignore_errors, prediction_response_df = combined)
  output$setting = setting$name
  output$EV_cargo = EV_cargo_oi
  # output$prediction_response_df = combined
  return(output)
}

#' @title Convert gene list to correct settings format for evaluation of target gene prediction.
#'
#' @description \code{convert_gene_list_settings_evaluation} Converts a gene list to correct settings format for evaluation of target gene prediction.
#'
#' @usage
#' convert_gene_list_settings_evaluation(gene_list, name, EV_cargo_oi, background)
#'
#' @param gene_list A character vector of target gene names
#' @param name The name that will be given to the setting
#' @param EV_cargo_oi The possibly active EV_cargo
#' @param background A character vector of names of genes that are not target genes. If genes present in the gene list are in this vector of background gene names, these genes will be removed from the background.

#' @return A list with following elements: $name, $from, $response
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' all_genes = unique(c(weighted_networks$gr$from,weighted_networks$gr$to,weighted_networks$lr_sig$from, weighted_networks$lr_sig$to))
#' gene_list = c("ID1","ID2","ID3")
#' setting = list(convert_gene_list_settings_evaluation(gene_list = c("ID1","ID2","ID3"), name = "test",EV_cargo_oi = "TGFB1", background = all_genes))
#' }
#' @export
#'
#'
convert_gene_list_settings_evaluation = function(gene_list, name, EV_cargo_oi, background) {
  # input check
  if(!is.character(gene_list))
    stop("gene_list should be character vector")
  if(!is.character(name) | length(name) > 1)
    stop("name should be character vector of length 1")
  if(!is.character(EV_cargo_oi))
    stop("EV_cargo_oi should be character vector")
  if(!is.character(background))
    stop("background should be character vector")

  requireNamespace("dplyr")

  background = background[(background %in% gene_list) == FALSE]

  background_logical = rep(FALSE,times = length(background))
  names(background_logical) = background
  gene_list_logical = rep(TRUE,times = length(gene_list))
  names(gene_list_logical) = gene_list
  response = c(background_logical,gene_list_logical)

  return(list(name = name, from = EV_cargo_oi, response = response))
}



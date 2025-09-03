#' @title Convert settings to correct settings format for EV_cargo prediction.
#'
#' @description \code{convert_settings_EV_cargo_prediction} Converts settings to correct settings format for EV_cargo activity prediction. In this prediction problem, EV_cargo (out of a set of possibly active EV_cargo) will be ranked based on feature importance scores. The format can be made suited for: 1) validation of EV_cargo activity state prediction by calculating individual feature importane scores or 2) feature importance based on models with embedded feature importance determination; applications in which EV_cargo need to be scores based on their possible upstream activity: 3) by calculating individual feature importane scores or 4) feature importance based on models with embedded feature importance determination.
#'
#' @usage
#' convert_settings_EV_cargo_prediction(settings, all_EV_cargo, validation = TRUE, single = TRUE)
#'
#' @param settings A list of lists. Eeach sublist contains the following elements: .$name: name of the setting; .$from: name(s) of the EV_cargo(s) active in the setting of interest; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest.
#' @param all_EV_cargo A character vector of possible EV_cargo that will be considered for the EV_cargo activity state prediction.
#' @param validation TRUE if seetings need to be prepared for validation of EV_cargo activity state predictions (this implies that the true active EV_cargo of a setting is known); FALSE for application purposes when the true active EV_cargo(s) is/are not known.
#' @param single TRUE if feature importance scores for EV_cargo will be calculated by looking at ligans individually. FALSE if the goal is to calculate the feature importance scores via sophisticated classification algorithms like random forest.

#' @return A list with following elements: $name, $EV_cargo: name of active EV_cargo(s) (only if validation is TRUE), $from (EV_cargo(s) that will be tested for activity prediction), $response
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
#' EV_cargo = unlist(extract_EV_cargo_from_settings(settings,combination = FALSE))
#' settings_EV_cargo_pred = convert_settings_EV_cargo_prediction(settings, EV_cargo, validation = TRUE, single = TRUE)
#' }
#' @export
#'
#'
convert_settings_EV_cargo_prediction = function(settings,all_EV_cargo,validation = TRUE, single = TRUE){

  # input check
  if(!is.list(settings))
    stop("settings should be a list")
  if(!is.character(all_EV_cargo))
    stop("all_EV_cargo should be a character vector")
  if(!is.logical(validation) | length(validation) != 1)
    stop("validation should be TRUE or FALSE")
  if(!is.logical(single) | length(single) != 1)
    stop("single should be TRUE or FALSE")

  requireNamespace("dplyr")

  new_settings = list()
  if (validation == TRUE && single == TRUE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      for (k in 1:length(all_EV_cargo)){
        test_EV_cargo = all_EV_cargo[[k]]
        new_settings[[length(new_settings) + 1]] = list(make_new_setting_EV_cargo_prediction_single_validation(setting,test_EV_cargo))
      }
    }
  } else if (validation == TRUE && single == FALSE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      new_settings[[length(new_settings) + 1]] = list(make_new_setting_EV_cargo_prediction_multi_validation(setting,all_EV_cargo))
    }
  } else if (validation == FALSE && single == TRUE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      for (k in 1:length(all_EV_cargo)){
        test_EV_cargo = all_EV_cargo[[k]]
        new_settings[[length(new_settings) + 1]] = list(make_new_setting_EV_cargo_prediction_single_application(setting,test_EV_cargo))
      }
    }
  } else if (validation == FALSE && single == FALSE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      new_settings[[length(new_settings) + 1]] = list(make_new_setting_EV_cargo_prediction_multi_application(setting,all_EV_cargo))
    }
  }
  return(new_settings %>% unlist(recursive = FALSE))
}
#' @title Get EV_cargo importances based on target gene prediction performance of single EV_cargo.
#'
#' @description \code{get_single_EV_cargo_importances} Get EV_cargo importance measures for EV_cargo based on how well a single, individual, EV_cargo can predict an observed response. Assess how well every EV_cargo of interest is able to predict the observed transcriptional response in a particular dataset, according to the EV_cargo-target model. It can be assumed that the EV_cargo that best predicts the observed response, is more likely to be the true EV_cargo.
#'
#' @usage
#' get_single_EV_cargo_importances(setting,EV_cargo_target_matrix, EV_cargo_position = "cols", known = TRUE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the EV_cargo(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $EV_cargo: NULL or the name of the EV_cargo(s) that are known to be active in the setting of interest.
#' @param known Indicate whether the true active EV_cargo for a particular dataset is known or not. Default: TRUE. The true EV_cargo will be extracted from the $EV_cargo slot of the setting.
#' @inheritParams evaluate_target_prediction
#'
#' @return A data.frame with for each EV_cargo - data set combination, classification evaluation metrics indicating how well the query EV_cargo predicts the response in the particular dataset. Evaluation metrics are the same as in \code{\link{evaluate_target_prediction}}. In addition to the metrics, the name of the particular setting ($setting), the name of the query EV_cargo($test_EV_cargo), the name of the true active EV_cargo (if known: $EV_cargo).
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_EV_cargo_pred = convert_settings_EV_cargo_prediction(settings, all_EV_cargo = unlist(extract_EV_cargo_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' EV_cargo = extract_EV_cargo_from_settings(settings_EV_cargo_pred,combination = FALSE)
#' EV_cargo_target_matrix = construct_EV_cargo_target_matrix(weighted_networks, EV_cargo)
#' EV_cargo_importances = dplyr::bind_rows(lapply(settings_EV_cargo_pred,get_single_EV_cargo_importances,EV_cargo_target_matrix))
#' print(head(EV_cargo_importances))
#' }
#' @export
#'
get_single_EV_cargo_importances = function(setting,EV_cargo_target_matrix, EV_cargo_position = "cols", known = TRUE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

  metrics = evaluate_target_prediction(setting, EV_cargo_target_matrix, EV_cargo_position)
  metrics = metrics %>% rename(test_EV_cargo = EV_cargo)
  if (known == TRUE){
    true_EV_cargo = setting$EV_cargo
    metrics_meta = metrics %>% select(setting,test_EV_cargo) %>% bind_cols(tibble(EV_cargo = true_EV_cargo))
    metrics = inner_join(metrics_meta, metrics, by = c("setting","test_EV_cargo"))
  }
  return(metrics)
}
#' @title Get EV_cargo importances from a multi-EV_cargo classfication model.
#'
#' @description \code{get_multi_EV_cargo_importances} A classificiation algorithm chosen by the user is trained to construct one model based on the target gene predictions of all EV_cargo of interest (EV_cargo are considered as features) in order to predict the observed response in a particular dataset. Variable importance scores that indicate for each EV_cargo the importance for response prediction, are extracted. It can be assumed that EV_cargo with higher variable importance scores are more likely to be a true active EV_cargo.
#'
#' @usage
#' get_multi_EV_cargo_importances(setting,EV_cargo_target_matrix, EV_cargo_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE, known = TRUE, filter_genes = FALSE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the EV_cargo(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $EV_cargo: NULL or the name of the EV_cargo(s) that are known to be active in the setting of interest.
#' @param EV_cargo_target_matrix A matrix of EV_cargo-target probabilty scores (recommended) or discrete target assignments (not-recommended).
#' @param EV_cargo_position Indicate whether the EV_cargo in the EV_cargo-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#' @param algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost"; without: "lda","naive_bayes","pls"(because bug in current version of pls package), "pcaNNet". Please notice that not all these algorithms work when the features (i.e. EV_cargo vectors) are categorical (i.e. discrete class assignments).
#' @param cv Indicate whether model training and hyperparameter optimization should be done via cross-validation. Default: TRUE. FALSE might be useful for applications only requiring variable importance, or when final model is not expected to be extremely overfit.
#' @param cv_number The number of folds for the cross-validation scheme: Default: 4; only relevant when cv == TRUE.
#' @param cv_repeats The number of repeats during cross-validation. Default: 2; only relevant when cv == TRUE.
#' @param parallel Indiciate whether the model training will occur parallelized. Default: FALSE. TRUE only possible for non-windows OS.
#' @param n_cores The number of cores used for parallelized model training via cross-validation. Default: 4. Only relevant on non-windows OS.
#' @param ignore_errors Indiciate whether errors during model training by caret should be ignored such that another model training try will be initiated until model is trained without raising errors. Default: FALSE.
#' @param continuous Indicate whether during training of the model, model training and evaluation should be done on class probabilities or discrete class labels. For huge class imbalance, we recommend setting this value to TRUE. Default: TRUE.
#' @param known Indicate whether the true active EV_cargo for a particular dataset is known or not. Default: TRUE. The true EV_cargo will be extracted from the $EV_cargo slot of the setting.
#' @param filter_genes Indicate whether 50 per cent of the genes that are the least variable in EV_cargo-target scores should be removed in order to reduce the training of the model. Default: FALSE.
#'
#' @return A data.frame with for each EV_cargo - data set combination, feature importance scores indicating how important the query EV_cargo is for the prediction of the response in the particular dataset, when prediction is done via a trained classification model with all possible EV_cargo as input. In addition to the importance score(s), the name of the particular setting ($setting), the name of the query EV_cargo($test_EV_cargo), the name of the true active EV_cargo (if known: $EV_cargo).
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_EV_cargo_pred = convert_settings_EV_cargo_prediction(settings, all_EV_cargo = unlist(extract_EV_cargo_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' EV_cargo = extract_EV_cargo_from_settings(settings_EV_cargo_pred,combination = FALSE)
#' EV_cargo_target_matrix = construct_EV_cargo_target_matrix(weighted_networks, EV_cargo)
#' EV_cargo_importances_glm = dplyr::bind_rows(lapply(settings_EV_cargo_pred, get_multi_EV_cargo_importances,EV_cargo_target_matrix, algorithm = "glm"))
#' print(head(EV_cargo_importances_glm))
#' }
#' @export
#'
get_multi_EV_cargo_importances = function(setting,EV_cargo_target_matrix, EV_cargo_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE, known = TRUE, filter_genes = FALSE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")
  if(!is.logical(filter_genes) | length(filter_genes) > 1)
    stop("filter_genes should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

  if (filter_genes == TRUE){
    EV_cargo_target_matrix = filter_genes_EV_cargo_target_matrix(EV_cargo_target_matrix,EV_cargo_position)
  }

  setting_name = setting$name
  output = evaluate_multi_EV_cargo_target_prediction(setting, EV_cargo_target_matrix, EV_cargo_position,algorithm, var_imps = TRUE, cv, cv_number, cv_repeats, parallel, n_cores, ignore_errors, continuous)
  metrics = output$var_imps
  metrics = metrics %>% mutate(setting = setting_name) %>% rename(test_EV_cargo = feature)

  if (known == TRUE){
    true_EV_cargo = setting$EV_cargo
    metrics = metrics %>% mutate(EV_cargo = true_EV_cargo)
    metrics = metrics %>% select(setting, test_EV_cargo, EV_cargo, importance)
    return(metrics)
  }
  metrics = metrics %>% select(setting, test_EV_cargo, importance)
  return(metrics)

}


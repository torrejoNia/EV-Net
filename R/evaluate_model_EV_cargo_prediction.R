#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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


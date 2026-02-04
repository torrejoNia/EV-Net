#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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



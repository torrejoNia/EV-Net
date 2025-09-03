
evaluate_target_prediction_strict = function(response,prediction,continuous = TRUE, prediction_response_df = FALSE){
  response_df = tibble(gene = names(response), response = response)
  prediction_df = tibble(gene = names(prediction), prediction = prediction)
  combined = inner_join(response_df,prediction_df, by = "gene")
  if (nrow(combined) == 0)
    stop("Gene names in response don't accord to gene names in EV_cargo-target matrix (did you consider differences human-mouse namings?)")
  prediction_vector = combined$prediction
  names(prediction_vector) = combined$gene
  response_vector = combined$response
  names(response_vector) = combined$gene
  if (continuous == TRUE){
    performance = classification_evaluation_continuous_pred(prediction_vector,response_vector)

  } else{
    performance = classification_evaluation_categorical_pred(prediction_vector,response_vector)
  }
  if (prediction_response_df == TRUE){
    output = list(performance = performance, prediction_response_df = combined)
    return(output)
  } else {
    return(performance)
  }

}

classification_evaluation_continuous_pred = function(prediction,response, iregulon = TRUE){

  if ((sd(response) == 0 & sd(prediction) == 0) | is.null(prediction) | is.null(response)){ # problems can occur otherwise in these leave-one-in models
    return(dplyr::tibble(auroc = NA,
                         aupr = NA,
                         aupr_corrected = NA,
                         sensitivity_roc = NA,
                         specificity_roc = NA,
                         mean_rank_GST_log_pval = NA,
                         auc_iregulon = NA,
                         auc_iregulon_corrected = NA,
                         pearson = NA,
                         spearman = NA))
  }
  prediction_ROCR = ROCR::prediction(prediction, response)
  performance1 = ROCR::performance(prediction_ROCR, measure="tpr", x.measure="fpr")

  performance2 = ROCR::performance(prediction_ROCR, measure="prec", x.measure="rec")
  performance = tibble(tpr = performance1@y.values[[1]], fpr=performance1@x.values[[1]], precision=performance2@y.values[[1]], recall=performance2@x.values[[1]])

  performance = performance %>% replace_na(list(recall=0, precision=1))

  aupr = caTools::trapz(performance$recall, performance$precision)
  pos_class = sum(response)
  total = length(response)
  aupr_random = pos_class/total


  metrics = get_split_auroc(prediction, response)

  sensitivity = metrics$auroc_sensitivity
  specificity = metrics$auroc_specificity
  auroc = metrics$auroc

  cor_p = cor(prediction, response)
  cor_s = cor(prediction, response, method = "s")

  cor_p_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response))) %>% .$p.value
  cor_s_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response), method =  "s")) %>% .$p.value

  # Mean rank GST calculated if limma is installed
  mean_rank_GST = ifelse(rlang::is_installed("limma"), limma::wilcoxGST(response, prediction), NA)

  # Calculate the AUC-iRegulon
  output_iregulon = list()
  if (iregulon){
    output_iregulon = calculate_auc_iregulon(prediction,response)
  }

  tbl_perf = tibble(auroc = auroc,
                    aupr = aupr,
                    aupr_corrected = aupr - aupr_random,
                    sensitivity_roc = sensitivity,
                    specificity_roc = specificity,
                    mean_rank_GST_log_pval = -log(mean_rank_GST),
                    auc_iregulon = output_iregulon$auc_iregulon,
                    auc_iregulon_corrected = output_iregulon$auc_iregulon_corrected,
                    pearson_log_pval = -log10(cor_p_pval),
                    spearman_log_pval = -log10(cor_s_pval),
                    pearson = cor_p,
                    spearman = cor_s)

  # Remove mean_rank_GST if limma is not installed
  if (!rlang::is_installed("limma")) {
    tbl_perf = tbl_perf %>% select(-mean_rank_GST_log_pval)
  }

  return(tbl_perf)
}

classification_evaluation_categorical_pred = function(predictions, response) {
  # print(head(predictions))
  # print(length(predictions))

  if (sd(response) == 0){ # if all response is the same
    return(dplyr::tibble(accuracy = NA,
                         recall = NA,
                         specificity = NA,
                         precision = NA,
                         F1 =  NA,
                         F05 = NA,
                         F2 = NA,
                         mcc = NA,
                         informedness = NA,
                         markedness = NA,
                         fisher_pval_log = NA,
                         fisher_odds = NA))
  }

  num_positives = sum(response)
  num_total = length(response)

  # calculate base statistics
  pos_preds = sum(predictions)

  tp = sum(response[predictions])
  fp = pos_preds - tp

  num_negatives = num_total - num_positives

  tp = tp
  fp = fp
  fn = num_positives - tp
  tn = num_negatives - fp
  npv = tn / (tn + fn)
  if (sd(predictions) == 0){
    fisher = list(p.value = NA, estimate = NA)
  } else {
    fisher = fisher.test(as.factor(response), predictions)
  }

  mcc_S = (tp + fn)/num_total
  mcc_P = (tp + fp)/num_total

  metrics = dplyr::tibble(
    accuracy = (tp + tn) / (num_positives + num_negatives),
    recall = tp / num_positives,
    specificity = tn / num_negatives,
    precision = tp / (tp + fp),
    F1 =  (2 * precision * recall) / (precision + recall),
    F05 = (1.25 * precision * recall) / (0.25 * precision + recall),
    F2 = (5 * precision * recall) / (4 * precision + recall),
    mcc = (tp/num_total - mcc_S * mcc_P)/sqrt(mcc_P * mcc_S * (1-mcc_S) * (1-mcc_P)) ,
    informedness = recall + specificity - 1,
    markedness = precision + npv - 1,
    fisher_pval_log = -log(fisher$p.value),
    fisher_odds = fisher$estimate
  )
  if (sd(predictions) == 0){ # all predictions are the same!
    return(dplyr::tibble(accuracy = metrics$accuracy,
                         recall = metrics$recall,
                         specificity = metrics$specificity,
                         precision = 0,
                         F1 =  0,
                         F05 = 0,
                         F2 = 0,
                         mcc = 0,
                         informedness = 0,
                         markedness = 0,
                         fisher_pval_log = NA,
                         fisher_odds = NA))
  }
  return(metrics)
}

make_new_setting_EV_cargo_prediction_single_validation = function(setting,test_EV_cargo){

  if (length(setting$from) > 1) {
    setting$from = paste0(setting$from,collapse = "-")
  }

  new_setting = list()
  new_setting$name = setting$name
  new_setting$EV_cargo = setting$from
  new_setting$from = test_EV_cargo
  new_setting$response = setting$response
  return(new_setting)
}

filter_genes_EV_cargo_target_matrix = function(EV_cargo_target_matrix, EV_cargo_position = cols){
  if (EV_cargo_position == "cols"){
    target_genes = rownames(EV_cargo_target_matrix)
    sd_genes = apply(EV_cargo_target_matrix,1,sd)
    EV_cargo_target_matrix_ = EV_cargo_target_matrix[sd_genes > quantile(sd_genes,0.5),]
  } else if (EV_cargo_position == "rows") {
    target_genes = colnames(EV_cargo_target_matrix)
    sd_genes = apply(EV_cargo_target_matrix,2,sd)
    EV_cargo_target_matrix_ = EV_cargo_target_matrix[,sd_genes > quantile(sd_genes,0.5)]
  }
  return(EV_cargo_target_matrix_)
}

construct_EV_cargo_signaling_df = function(EV_cargo_all,targets_all,k,weighted_networks,EV_cargo_tf_matrix){
  final_combined_df = bind_rows(expand.grid(EV_cargo_all,targets_all) %>% apply(.,1,wrappper_visualization,k,weighted_networks,EV_cargo_tf_matrix))
}

wrappper_visualization = function(grid,k,weighted_networks,EV_cargo_tf_matrix){
  EV_cargo = grid[1]
  target = grid[2]
  network = get_network_df(EV_cargo,target,k,weighted_networks,EV_cargo_tf_matrix)
}

get_network_df = function(EV_cargo_to_vis,target_to_vis,k,weighted_networks,EV_cargo_tf_matrix){

  ## prepare TFs downstream of EV_cargo

  ## EV_cargo should be in columns

  EV_cargo_tf_matrix_visualization = EV_cargo_tf_matrix[,EV_cargo_to_vis] %>% as.matrix()

  colnames(EV_cargo_tf_matrix_visualization) = "V1"

  EV_cargo_tf_matrix_visualization_df = as_tibble(EV_cargo_tf_matrix_visualization) %>% rename(weight = V1) %>% mutate(TF = rownames(EV_cargo_tf_matrix)) %>% select(TF,weight)
  EV_cargo_tf_matrix_visualization_df_filtered = EV_cargo_tf_matrix_visualization_df %>% filter(weight > 0) %>% mutate(EV_cargo = EV_cargo_to_vis)

  ## prepare TFs upstream of target
  regulatory_network_filtered = weighted_networks$gr %>% filter(to == target_to_vis) %>% rename(TF = from, weight_grn = weight)

  ## combine both
  combined_df = inner_join(EV_cargo_tf_matrix_visualization_df_filtered,regulatory_network_filtered, by = "TF")
  combined_df = combined_df %>% mutate(total_weight = weight*weight_grn)
  final_combined_df = combined_df %>% arrange(-total_weight) %>% .[1:min(k,nrow(combined_df)),]

  return(final_combined_df)
}

make_new_setting_EV_cargo_prediction_single_application = function(setting,test_EV_cargo){
  new_setting = list()
  new_setting$name = setting$name
  new_setting$from = test_EV_cargo
  new_setting$response = setting$response
  return(new_setting)
}

get_split_auroc = function(observed, known) {
  prediction = ROCR::prediction(observed, known)
  performance = ROCR::performance(prediction, measure="spec", x.measure="rec")
  metrics = tibble(tpr = performance@x.values[[1]], spec = performance@y.values[[1]])
  meetingpoint = which(-(1-metrics$spec) + 1 <= metrics$tpr)[[1]] # < or <= ?
  specs = c((1-metrics$spec)[seq_len(meetingpoint)-1],1-metrics$tpr[meetingpoint], 1)
  recs = c(metrics$tpr[seq_len(meetingpoint)], 0)
  auroc_specificity = caTools::trapz(specs, recs)
  auroc = caTools::trapz(1-metrics$spec, metrics$tpr)
  auroc_sensitivity = auroc - auroc_specificity
  tibble(auroc=auroc, auroc_sensitivity=auroc_sensitivity, auroc_specificity=auroc_specificity)
}

calculate_auc_iregulon = function(prior,response){

  genes_prior = names(prior)
  dim(prior) = c(1,length(prior))
  colnames(prior) = genes_prior
  rownames(prior) = "EV_cargo"

  prior_rank = apply(prior,1,rank_desc)
  rankings = tibble(prior=prior_rank[,1], rn = rownames(prior_rank))

  fake_rankings = rankings %>% mutate(rn = sample(rn))
  rankings = data.table::data.table(rankings)
  fake_rankings = data.table::data.table(fake_rankings)

  aucMaxRank = 0.03*nrow(rankings)

  # calculate enrichment over the expression settings

  geneSet = response[response == TRUE] %>% names()
  geneSet = unique(geneSet)
  nGenes = length(geneSet)
  geneSet = geneSet[which(geneSet %in% rankings$rn)]

  missing = nGenes-length(geneSet)

  gSetRanks = subset(rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed

  aucThreshold = round(aucMaxRank)
  maxAUC = aucThreshold * nrow(gSetRanks)

  auc_iregulon = sapply(gSetRanks, .calcAUC, aucThreshold, maxAUC)

  gSetRanks_fake = subset(fake_rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed
  auc_iregulon_fake = sapply(gSetRanks_fake, .calcAUC, aucThreshold, maxAUC)

  auc_iregulon_corrected = auc_iregulon - auc_iregulon_fake
  return(list(auc_iregulon = auc_iregulon, auc_iregulon_corrected = auc_iregulon_corrected))
}

rank_desc = function(x){rank(desc(x), ties.method = "max")}
# rank_desc = function(x){rank(desc(x), ties.method = "random")}

.calcAUC = function(oneRanking, aucThreshold, maxAUC)
{
  x = unlist(oneRanking)
  x = sort(x[x<aucThreshold])
  y = 1:length(x)
  a = diff(c(x, aucThreshold)) * y
  return(sum(a)/maxAUC)
}

get_shortest_path_signaling = function(EV_cargo_oi, signaling_df, signaling_igraph){
  EV_cargo_signaling = signaling_df %>% filter(EV_cargo == EV_cargo_oi)
  tfs = EV_cargo_signaling$TF %>% unique()
  sp = igraph::shortest_paths(signaling_igraph, from = EV_cargo_oi, to = tfs, mode ="out")
  tf_nodes = unlist(sp$vpath) %>% names() %>% unique() %>% .[. != EV_cargo_oi] # EV_cargo should not belong to tf nodes
}

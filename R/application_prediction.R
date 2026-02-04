#' @title Predict activities of EV_cargo in regulating expression of a gene set of interest
#'
#' @description \code{predict_EV_cargo_activities} Predict activities of EV_cargo in regulating expression of a gene set of interest. Ligand activities are defined as how well they predict the observed transcriptional response (i.e. gene set) according to the NicheNet model.
#'
#' @usage
#' predict_EV_cargo_activities(geneset, background_expressed_genes,EV_cargo_target_matrix, potential_EV_cargo, single = TRUE,...)
#'
#' @param geneset Character vector of the gene symbols of genes of which the expression is potentially affected by EV_cargo from the interacting cell.
#' @param background_expressed_genes Character vector of gene symbols of the background, non-affected, genes (can contain the symbols of the affected genes as well).
#' @param EV_cargo_target_matrix The NicheNet EV_cargo-target matrix denoting regulatory potential scores between EV_cargo and targets (EV_cargo in columns).
#' @param potential_EV_cargo Character vector giving the gene symbols of the potentially active EV_cargo you want to define EV_cargo activities for.
#' @param single TRUE if you want to calculate EV_cargo activity scores by considering every EV_cargo individually (recommended). FALSE if you want to calculate EV_cargo activity scores as variable importances of a multi-EV_cargo classification model.
#' @param ... Additional parameters for get_multi_EV_cargo_importances if single = FALSE.
#'
#' @return A tibble giving several EV_cargo activity scores. Following columns in the tibble: $test_EV_cargo, $auroc, $aupr and $pearson.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' EV_cargo = list("TNF","BMP2","IL4")
#' EV_cargo_target_matrix = construct_EV_cargo_target_matrix(weighted_networks, EV_cargo, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_EV_cargo = c("TNF","BMP2","IL4")
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' EV_cargo_activities = predict_EV_cargo_activities(geneset = geneset, background_expressed_genes = background_expressed_genes, EV_cargo_target_matrix = EV_cargo_target_matrix, potential_EV_cargo = potential_EV_cargo)
#' }
#'
#' @export
#'
predict_EV_cargo_activities = function(geneset,background_expressed_genes,EV_cargo_target_matrix, potential_EV_cargo, single = TRUE,...){
  setting = list(geneset) %>%
    lapply(convert_gene_list_settings_evaluation, name = "gene set", EV_cargo_oi = potential_EV_cargo, background = background_expressed_genes)
  if (single == TRUE){
    settings_EV_cargo_prediction = setting %>%
      convert_settings_EV_cargo_prediction(all_EV_cargo = potential_EV_cargo, validation = FALSE, single = TRUE)
    EV_cargo_importances = settings_EV_cargo_prediction %>% lapply(get_single_EV_cargo_importances,EV_cargo_target_matrix = EV_cargo_target_matrix, known = FALSE) %>% bind_rows()

  } else {
    settings_EV_cargo_prediction = setting %>%
      convert_settings_EV_cargo_prediction(all_EV_cargo = potential_EV_cargo, validation = FALSE, single = FALSE)
    EV_cargo_importances = settings_EV_cargo_prediction %>% lapply(get_multi_EV_cargo_importances,EV_cargo_target_matrix = EV_cargo_target_matrix, known = FALSE, ...) %>% bind_rows()

  }
  return(EV_cargo_importances %>% select(test_EV_cargo,auroc,aupr,aupr_corrected, pearson))
}
#' @title Infer weighted active EV_cargo-target links between a possible EV_cargo and target genes of interest
#'
#' @description \code{get_weighted_EV_cargo_target_links} Infer active EV_cargo target links between possible EV_cargo and genes belonging to a gene set of interest: consider the intersect between the top n targets of a EV_cargo and the gene set.
#'
#' @usage
#' get_weighted_EV_cargo_target_links(EV_cargo, geneset,EV_cargo_target_matrix,n = 250)
#'
#' @param geneset Character vector of the gene symbols of genes of which the expression is potentially affected by EV_cargo from the interacting cell.
#' @param EV_cargo Character vector giving the gene symbols of the potentially active EV_cargo for which you want to find target genes.
#' @param n The top n of targets per EV_cargo that will be considered. Default: 250.
#' @inheritParams predict_EV_cargo_activities
#'
#' @return A tibble with columns EV_cargo, target and weight (i.e. regulatory potential score).
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' EV_cargo = list("TNF","BMP2","IL4")
#' EV_cargo_target_matrix = construct_EV_cargo_target_matrix(weighted_networks, EV_cargo, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_EV_cargo = "TNF"
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' active_EV_cargo_target_links_df = get_weighted_EV_cargo_target_links(EV_cargo = potential_EV_cargo, geneset = geneset, EV_cargo_target_matrix = EV_cargo_target_matrix, n = 250)
#' }
#'
#' @export
#'
get_weighted_EV_cargo_target_links = function(EV_cargo, geneset,EV_cargo_target_matrix,n = 250){
  top_n_score = EV_cargo_target_matrix[,EV_cargo] %>% sort(decreasing = T) %>% head(n) %>% min()
  targets = intersect(EV_cargo_target_matrix[,EV_cargo] %>% .[. >= top_n_score ] %>% names(),geneset)
  if (length(targets) == 0){
    EV_cargo_target_weighted_df = tibble(EV_cargo = EV_cargo, target = NA, weight = NA)
  } else if (length(targets) == 1) {
    EV_cargo_target_weighted_df = tibble(EV_cargo = EV_cargo, target = targets, weight = EV_cargo_target_matrix[targets,EV_cargo])
  } else {
    EV_cargo_target_weighted_df = tibble(EV_cargo = EV_cargo, target = names(EV_cargo_target_matrix[targets,EV_cargo])) %>% inner_join(tibble(target = names(EV_cargo_target_matrix[targets,EV_cargo]), weight = EV_cargo_target_matrix[targets,EV_cargo]), by = "target")
  }
  return(EV_cargo_target_weighted_df)
}
#' @title Prepare heatmap visualization of the EV_cargo-target links starting from a EV_cargo-target tibble.
#'
#' @description \code{prepare_EV_cargo_target_visualization} Prepare heatmap visualization of the EV_cargo-target links starting from a EV_cargo-target tibble. Get regulatory potential scores between all pairs of EV_cargo and targets documented in this tibble. For better visualization, we propose to define a quantile cutoff on the EV_cargo-target scores.
#'
#' @usage
#' prepare_EV_cargo_target_visualization(EV_cargo_target_df, EV_cargo_target_matrix, cutoff = 0.25)
#'
#' @param cutoff Quantile cutoff on the EV_cargo-target scores of the input weighted EV_cargo-target network. Scores under this cutoff will be set to 0.
#' @param EV_cargo_target_df Tibble with columns 'EV_cargo', 'target' and 'weight' to indicate EV_cargo-target regulatory potential scores of interest.
#' @inheritParams predict_EV_cargo_activities
#'
#' @return A matrix giving the EV_cargo-target regulatory potential scores between EV_cargo of interest and their targets genes part of the gene set of interest.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' EV_cargo = list("TNF","BMP2","IL4")
#' EV_cargo_target_matrix = construct_EV_cargo_target_matrix(weighted_networks, EV_cargo, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' active_EV_cargo_target_links_df = potential_EV_cargo %>% lapply(get_weighted_EV_cargo_target_links, geneset = geneset, EV_cargo_target_matrix = EV_cargo_target_matrix, n = 250) %>% bind_rows()
#' active_EV_cargo_target_links = prepare_EV_cargo_target_visualization(EV_cargo_target_df = active_EV_cargo_target_links_df, EV_cargo_target_matrix = EV_cargo_target_matrix, cutoff = 0.25)
#' }
#'
#' @export
#'
prepare_EV_cargo_target_visualization = function(EV_cargo_target_df, EV_cargo_target_matrix, cutoff = 0.25){

  # define a cutoff on the EV_cargo-target links
  cutoff_include_all_EV_cargo = EV_cargo_target_df$weight %>% quantile(cutoff)

  # give a score of 0 to EV_cargo-target links not higher than the defined cutoff
  EV_cargo_target_matrix_oi = EV_cargo_target_matrix
  EV_cargo_target_matrix_oi[EV_cargo_target_matrix_oi < cutoff_include_all_EV_cargo] = 0

  # consider only targets belonging to the top250 targets of individual EV_cargo and with at least one EV_cargo-link with score higher than the defined cutoff
  EV_cargo_target_vis = EV_cargo_target_matrix_oi[EV_cargo_target_df$target %>% unique(),EV_cargo_target_df$EV_cargo %>% unique()]
  dim(EV_cargo_target_vis) = c(length(EV_cargo_target_df$target %>% unique()), length(EV_cargo_target_df$EV_cargo %>% unique()))
  all_targets = EV_cargo_target_df$target %>% unique()
  all_EV_cargo = EV_cargo_target_df$EV_cargo %>% unique()
  rownames(EV_cargo_target_vis) = all_targets
  colnames(EV_cargo_target_vis) = all_EV_cargo

  keep_targets = all_targets[EV_cargo_target_vis %>% apply(1,sum) > 0]
  keep_EV_cargo = all_EV_cargo[EV_cargo_target_vis %>% apply(2,sum) > 0]


  EV_cargo_target_vis_filtered = EV_cargo_target_vis[keep_targets,keep_EV_cargo]


  if(is.matrix(EV_cargo_target_vis_filtered)){
    rownames(EV_cargo_target_vis_filtered) = keep_targets
    colnames(EV_cargo_target_vis_filtered) = keep_EV_cargo

  } else {
    dim(EV_cargo_target_vis_filtered) = c(length(keep_targets), length(keep_EV_cargo))
    rownames(EV_cargo_target_vis_filtered) = keep_targets
    colnames(EV_cargo_target_vis_filtered) = keep_EV_cargo
  }

  if(nrow(EV_cargo_target_vis_filtered) > 1 & ncol(EV_cargo_target_vis_filtered) > 1){
    distoi = dist(1-cor(t(EV_cargo_target_vis_filtered)))
    hclust_obj = hclust(distoi, method = "ward.D2")
    order_targets = hclust_obj$labels[hclust_obj$order]

    distoi_targets = dist(1-cor(EV_cargo_target_vis_filtered))
    hclust_obj = hclust(distoi_targets, method = "ward.D2")
    order_EV_cargo = hclust_obj$labels[hclust_obj$order]

  } else {
    order_targets = rownames(EV_cargo_target_vis_filtered)
    order_EV_cargo = colnames(EV_cargo_target_vis_filtered)
  }

  vis_EV_cargo_target_network = EV_cargo_target_vis_filtered[order_targets,order_EV_cargo]
  dim(vis_EV_cargo_target_network) = c(length(order_targets), length(order_EV_cargo))
  rownames(vis_EV_cargo_target_network) = order_targets
  colnames(vis_EV_cargo_target_network) = order_EV_cargo
  return(vis_EV_cargo_target_network)

}

#' @keywords internal
#' 
get_weighted_EV_cargo_receptor_links = function(best_upstream_EV_cargo, expressed_receptors, lr_network, weighted_networks_lr_sig) {

  lr_network <- lr_network %>% distinct(from, to)
  weighted_networks_lr <- inner_join(weighted_networks_lr_sig, lr_network, by = c("from","to"))

  lr_network_top <- lr_network %>% filter(from %in% best_upstream_EV_cargo & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

  lr_network_top_df_long <- weighted_networks_lr %>% filter(from %in% best_upstream_EV_cargo & to %in% best_upstream_receptors)

  return(lr_network_top_df_long)

}


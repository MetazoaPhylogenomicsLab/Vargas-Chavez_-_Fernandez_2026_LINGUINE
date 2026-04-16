#' Create LINGUINE Configuration Object
#'
#' @description Generates the strict parameter list required to run the pipeline.
#'
#' @param dataset Character. Name of the dataset/run.
#' @param min_chromosome_length_bp Numeric. Minimum scaffold size.
#' @param orthology_type Character. "OGs" or "HOGs".
#' @param orthology_filename Character. Target file or directory name.
#' @param tree_filename Character. Newick tree file name.
#' @param base_dir Character. Root directory for the run.
#' @param raw_data_dir Character. Optional absolute path to raw data. Defaults to base_dir/raw_data.
#' @param prob_self_syn_hmm1 Numeric. HMM self-transition probability.
#' @param prob_to_other_syn_hmm1 Numeric. HMM distinct-transition probability.
#' @param prob_emit_correct_syn_hmm1_value Numeric. HMM correct emission probability.
#' @param prob_emit_non_syn_obs_from_syn_hmm1_value Numeric. HMM noise emission probability.
#' @param prob_emit_no_ortholog_from_syn_hmm1_value Numeric. HMM missing emission probability.
#' @param prob_emit_multiple_b_chrs_from_syn_hmm1_value Numeric. HMM paralogy emission probability.
#' @param window_size Numeric. Sliding window size for block classification.
#' @param purity_threshold Numeric. Minimum gene purity for block consensus.
#' @param min_final_block_size Numeric. Minimum genes to retain a block.
#' @param fisher_p Numeric. Enrichment p-value threshold.
#' @param fisher_odds Numeric. Enrichment odds ratio threshold.
#' @param min_ogs_by_block Numeric. Minimum distinct OGs to validate a block.
#' @param og_abundance Numeric. Percentile threshold to discard ubiquitous OGs.
#' @param fisher_p_lgs Numeric. Graph clustering p-value threshold.
#' @param fisher_odds_lgs Numeric. Graph clustering odds ratio threshold.
#' @param paralogy_p Numeric. Paralogy consolidation p-value.
#' @param paralogy_odds Numeric. Paralogy consolidation odds ratio.
#' @param outgroup_dominance Numeric. Minimum proportion for outgroup structural candidate.
#' @param parent_assignment Numeric. Minimum overlap for orphan partition assignment.
#' @param ribbon_max_links Numeric. Maximum visual links per orthogroup.
#' @param min_lg_fraction Numeric. Minimum fraction of total orthogroups an LG must contain to be retained (e.g., 0.01 for 1%).
#' @param resolve_multimapped Character. Strategy for handling orthogroups mapped to multiple LGs. Options: "drop" (strict 1-to-1), "random" (preserve density), "keep" (preserve paralogy).
#'
#' @return A `linguine_config` list object.
#' @export
create_linguine_config <- function(
    dataset,
    min_chromosome_length_bp = 0,
    orthology_type = "OGs",
    orthology_filename = "Orthogroups.tsv",
    tree_filename = "species_tree.nwk",
    base_dir = getwd(),
    raw_data_dir = NULL,
    prob_self_syn_hmm1 = 0.995,
    prob_to_other_syn_hmm1 = 0.005,
    prob_emit_correct_syn_hmm1_value = 0.85,
    prob_emit_non_syn_obs_from_syn_hmm1_value = 0.08,
    prob_emit_no_ortholog_from_syn_hmm1_value = 0.05,
    prob_emit_multiple_b_chrs_from_syn_hmm1_value = 0.01,
    window_size = 10,
    purity_threshold = 0.80,
    min_final_block_size = 10,
    fisher_p = 0.01,
    fisher_odds = 2.0,
    min_ogs_by_block = 2,
    og_abundance = 90,
    fisher_p_lgs = 1e-20,
    fisher_odds_lgs = 4,
    paralogy_p = 1e-20,
    paralogy_odds = 10.0,
    outgroup_dominance = 0.05,
    parent_assignment = 0.60,
    ribbon_max_links = 5,
    min_lg_fraction = 0.01,
    resolve_multimapped = "drop"
) {

  run_id <- paste0(dataset, "_min_chr_size_", format(min_chromosome_length_bp, scientific = FALSE), "bp")
  raw_dir_path <- if (!is.null(raw_data_dir)) raw_data_dir else file.path(base_dir, "raw_data")

  config <- list(
    dataset = dataset,
    run_id = run_id,
    min_chromosome_length_bp = min_chromosome_length_bp,
    orthology_type = orthology_type,
    orthology_filename = orthology_filename,
    tree_filename = tree_filename,

    paths = list(
      raw_data = raw_dir_path,
      processed_data = file.path(base_dir, "runs", run_id, "processed_data"),
      intermediate_data = file.path(base_dir, "runs", run_id, "intermediate_data"),
      results = file.path(base_dir, "runs", run_id, "results"),
      plots = file.path(base_dir, "runs", run_id, "plots")
    ),

    hmm = list(
      prob_self_syn_hmm1 = prob_self_syn_hmm1,
      prob_to_other_syn_hmm1 = prob_to_other_syn_hmm1,
      prob_emit_correct_syn_hmm1_value = prob_emit_correct_syn_hmm1_value,
      prob_emit_non_syn_obs_from_syn_hmm1_value = prob_emit_non_syn_obs_from_syn_hmm1_value,
      prob_emit_no_ortholog_from_syn_hmm1_value = prob_emit_no_ortholog_from_syn_hmm1_value,
      prob_emit_multiple_b_chrs_from_syn_hmm1_value = prob_emit_multiple_b_chrs_from_syn_hmm1_value
    ),

    thresholds = list(
      window_size = window_size,
      purity_threshold = purity_threshold,
      min_final_block_size = min_final_block_size,
      fisher_p = fisher_p,
      fisher_odds = fisher_odds,
      min_ogs_by_block = min_ogs_by_block,
      og_abundance = og_abundance,
      fisher_p_lgs = fisher_p_lgs,
      fisher_odds_lgs = fisher_odds_lgs,
      paralogy_p = paralogy_p,
      paralogy_odds = paralogy_odds,
      outgroup_dominance = outgroup_dominance,
      parent_assignment = parent_assignment,
      ribbon_max_links = ribbon_max_links,
      min_lg_fraction = min_lg_fraction,
      resolve_multimapped = resolve_multimapped
    )
  )
  class(config) <- "linguine_config"
  return(config)
}

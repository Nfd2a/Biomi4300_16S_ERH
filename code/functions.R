
matround <- function(x){trunc(x+0.5)}

scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "round") {

  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, function(x) {(n * x/sum(x))})

  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- ceiling(otu_table(physeq.scale))
  } else if (round == "matround"){
    otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
  } else if (round == "exactround"){
    # exactround calculates the rounding threshold for each sample so that 
    # all sample sums normalize exactly to the previous minimum sample sum
    transposed_otu <- t(otu_table(physeq.scale))
    for(sample in 1:nrow(transposed_otu)){
      sample_vec <- transposed_otu[sample,]
      asvs_to_round_up <- n - sum(floor(sample_vec))
      if (asvs_to_round_up == 0)
      {
        next
      }
      sample_sorted <- sort(sapply(sample_vec, function(x){x %% 1}))
      threshold <- sample_sorted[length(sample_sorted) - asvs_to_round_up]
      round_on_threshold_arr <- sample(which((sample_vec - floor(sample_vec)) == threshold), length(which(which(sample_sorted == threshold) > (length(sample_sorted) - asvs_to_round_up))), replace = FALSE)
      exactround <- Vectorize(function(x, ind){
        if ((x - floor(x)) > threshold)
          {ceiling(x)} 
        else if ((x - floor(x)) < threshold)
          {floor(x)}
        else {if (ind %in% round_on_threshold_arr)
            {
            return(ceiling(x))}
          else {
            return(floor(x))}}})
      transposed_otu[sample,] <- exactround(sample_vec, array(1:length(sample_vec)))
    }
    otu_table(physeq.scale) <- t(transposed_otu)
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}
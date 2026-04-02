validate_kmers <- function(kmers, k = 5) {
  if (!is.character(kmers)) stop("top_kmers must be a character vector.")
  if (length(kmers) < 1) stop("top_kmers is empty.")
  bad_len <- which(nchar(kmers) != k)
  if (length(bad_len) > 0) {
    stop(sprintf("All kmers must have length %d. Bad indices: %s",
                 k, paste(head(bad_len, 20), collapse = ",")))
  }
  if (any(!grepl("^[ACGT]+$", kmers))) {
    stop("All kmers must contain only A/C/G/T (uppercase).")
  }
  invisible(TRUE)
}


base_prob_from_kmers <- function(kmers, alpha = 1) {
  s <- paste0(kmers, collapse = "")
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  tab <- table(chars)
  
  bases <- c("A","C","G","T")
  counts <- setNames(rep(0, 4), bases)
  counts[names(tab)] <- as.numeric(tab)
  
  counts <- counts + alpha
  counts / sum(counts)
}


rand_dna <- function(n, prob = c(A=0.25, C=0.25, G=0.25, T=0.25)) {
  if (n <= 0) return("")
  paste0(sample(names(prob), size = n, replace = TRUE, prob = prob), collapse = "")
}

choose_overlap_soft <- function(contig, kmer, max_ol = 4,
                                mode = c("uniform", "biased", "score"),
                                beta = 0.6,
                                gamma = 1.5) {
  mode <- match.arg(mode)
  
  max_ol <- min(max_ol, nchar(kmer) - 1L, nchar(contig))
  if (max_ol <= 0) return(0L)
  
  ols <- 0:max_ol
  
  if (mode == "uniform") {
    return(sample(ols, 1))
  }
  
  if (mode == "biased") {

    w <- exp(beta * ols)
    w <- w / sum(w)
    return(sample(ols, 1, prob = w))
  }
  

  score <- vapply(ols, function(ol) {
    if (ol == 0) return(0)
    suf <- substr(contig, nchar(contig) - ol + 1, nchar(contig))
    pre <- substr(kmer, 1, ol)
    suf_c <- strsplit(suf, "", fixed = TRUE)[[1]]
    pre_c <- strsplit(pre, "", fixed = TRUE)[[1]]
    sum(suf_c == pre_c)
  }, numeric(1))
  

  w <- (score + 1) ^ gamma * exp(0.2 * ols)
  w <- w / sum(w)
  sample(ols, 1, prob = w)
}


assemble_from_4kmers <- function(top_kmers,
                                 target_len,
                                 k = 5,
                                 n_kmers = 4,
                                 max_ol = 4,
                                 overlap_mode = c("uniform", "biased", "score"),
                                 beta = 0.6,
                                 gamma = 1.5,
                                 start_random = TRUE,
                                 use_kmer_base_prob = TRUE,
                                 alpha = 1,
                                 max_attempts = 50,
                                 allow_repeat = FALSE) {
  
  overlap_mode <- match.arg(overlap_mode)

  stopifnot(is.numeric(target_len), length(target_len) == 1, target_len > 0)
  stopifnot(n_kmers >= 1, n_kmers <= length(top_kmers))
  stopifnot(max_ol >= 0, max_ol <= (k - 1))
  
  validate_kmers(top_kmers, k = k)
  
  for (attempt in seq_len(max_attempts)) {
    
    picked <- sample(top_kmers, size = n_kmers, replace = allow_repeat)
    if (start_random) picked <- sample(picked, length(picked))
    
    contig <- picked[1]
    overlaps_used <- integer(0)
    
    if (n_kmers > 1) {
      for (j in 2:n_kmers) {
        km <- picked[j]
        
        ol <- choose_overlap_soft(
          contig = contig,
          kmer = km,
          max_ol = max_ol,
          mode = overlap_mode,
          beta = beta,
          gamma = gamma
        )
        
        overlaps_used <- c(overlaps_used, ol)
        
   
        contig <- paste0(contig, substr(km, ol + 1, k))
      }
    }
    
    raw_len <- nchar(contig)
    

    if (raw_len > target_len) next
    
    fill_n <- target_len - raw_len
    if (fill_n > 0) {
      prob <- if (use_kmer_base_prob) {
        base_prob_from_kmers(picked, alpha = alpha)
      } else {
        c(A=0.25, C=0.25, G=0.25, T=0.25)
      }
      contig <- paste0(contig, rand_dna(fill_n, prob = prob))
    }
    
    return(list(
      contig_target = contig,
      raw_len = raw_len,
      fill_added = fill_n,
      picked_kmers = picked,
      overlaps = overlaps_used
    ))
  }
  

  picked <- sample(top_kmers, size = n_kmers, replace = allow_repeat)
  if (start_random) picked <- sample(picked, length(picked))
  contig <- paste0(picked, collapse = "")
  contig <- substr(contig, 1, target_len)
  
  list(
    contig_target = contig,
    raw_len = target_len,
    fill_added = 0L,
    picked_kmers = picked,
    overlaps = NA_integer_
  )
}


generate_large_batch_parallel <- function(top_kmers,
                                          lengths = 15:20,
                                          n_per_length = 200000,
                                          chunk_size = 20000,
                                          # old args kept for compatibility; no longer needed by the new assembler:
                                          sample_n = 4,
                                          min_ol = 1,
                                          random_tie = TRUE,
                                          start_random = TRUE,
                                          max_fill_random = 10,
                                          seed = 1,
                                          # new args:
                                          overlap_mode = c("uniform", "biased", "score"),
                                          max_ol = 4,
                                          beta = 0.6,
                                          gamma = 1.5,
                                          use_kmer_base_prob = TRUE,
                                          alpha = 1,
                                          max_attempts = 50,
                                          allow_repeat = FALSE,
                                          return_trace = FALSE) {
  
  overlap_mode <- match.arg(overlap_mode)
  
  validate_kmers(top_kmers, k = 5)
  
  task_df <- do.call(rbind, lapply(lengths, function(L) {
    starts <- seq(1, n_per_length, by = chunk_size)
    ends <- pmin(starts + chunk_size - 1, n_per_length)
    data.frame(target_len = L, start = starts, end = ends)
  }))
  task_df$task_id <- seq_len(nrow(task_df))
  
  blocks <- future.apply::future_lapply(task_df$task_id, function(tid) {
    
    L <- task_df$target_len[tid]
    s <- task_df$start[tid]
    e <- task_df$end[tid]
    n_block <- e - s + 1
    

    set.seed(seed + tid)
    
    seqs <- character(n_block)
    rawl <- integer(n_block)
    fill <- integer(n_block)
    
    if (return_trace) {
      picked_str <- character(n_block)
      ol_str <- character(n_block)
    }
    
    for (i in seq_len(n_block)) {
      out <- assemble_from_4kmers(
        top_kmers = top_kmers,
        target_len = L,
        k = 5,
        n_kmers = 4,
        max_ol = max_ol,
        overlap_mode = overlap_mode,
        beta = beta,
        gamma = gamma,
        start_random = start_random,
        use_kmer_base_prob = use_kmer_base_prob,
        alpha = alpha,
        max_attempts = max_attempts,
        allow_repeat = allow_repeat
      )
      
      seqs[i] <- out$contig_target
      rawl[i] <- out$raw_len
      fill[i] <- out$fill_added
      
      if (return_trace) {
        picked_str[i] <- paste(out$picked_kmers, collapse = ",")
        ol_str[i] <- paste(out$overlaps, collapse = ",")
      }
    }
    
    df <- data.frame(
      target_len = rep(L, n_block),
      sequence   = seqs,
      raw_len    = rawl,
      fill_added = fill,
      stringsAsFactors = FALSE
    )
    
    if (return_trace) {
      df$picked_kmers <- picked_str
      df$overlaps <- ol_str
    }
    
    df
  })
  
  big_df <- do.call(rbind, blocks)
  rownames(big_df) <- NULL
  big_df
}

data <- read.csv("combined_5mer_analysis_top80pct.csv")
data <- data[order(data$rank_mean,decreasing = F),]
kmers <- data$kmer
top_kmers <- kmers[1:100]
bot_kmers <- kmers[925:1024]

 library(future.apply)
 future::plan(future::multisession, workers = 63)
 df <- generate_large_batch_parallel(
   top_kmers,
   lengths = 15:20,
   n_per_length = 500000,
   chunk_size = 200,
   seed = 1,
   overlap_mode = "uniform", 
   return_trace = TRUE
 )

write.csv(df,"../Dataset/pos.csv")
df2 <- generate_large_batch_parallel(
   top_kmers,
   lengths = 15:20,
   n_per_length = 100000,
   chunk_size = 200,
   seed = 1,
   overlap_mode = "uniform", 
   return_trace = TRUE
)

write.csv(df2,"../Dataset/neg.csv")
 
 
 
#' Grammar Compressed Matrix Class
#'
#' @export
GrammarCompressedMatrix <- R6::R6Class(
  "GrammarCompressedMatrix",
  private = list(
    .ncol = 0L,
    .nrow = 0L,
    .compressed_matrix = NULL,
    .grammer_trees = NULL,
    .terminal_symbols = NULL,
    .P = NULL, # summation of terminal symbols
    .R = c(),
    .hash_table = character(),
    .to_diff_index_list = function(matrix) {
      index_list <- apply(matrix, 1, function(row) which(row == 1))
      diff_index_list <- lapply(index_list, function(vec) {
        if (length(vec) == 0) {
          integer(0)
        } else {
          c(vec[1], diff(vec))
        }
      })
    },
    .replace = function(diff_index_list, replace_pair_strs, new_symbols) {
      for (i in seq_along(replace_pair_strs)) {
        private$.hash_table[replace_pair_strs[i]] <- new_symbols[i]
      }
      replaces <- decode_pairs(replace_pair_strs)
      first_symbols <- vapply(replaces, function(pair) pair[1], integer(1))
      second_symbols <- vapply(replaces, function(pair) pair[2], integer(1))
      lapply(diff_index_list, function(vec) {
        n <- length(vec)
        if (n <= 1) return(vec)
        result <- integer(n)
        res_ind <- 1

        skip <- FALSE
        for (i in seq_len(n-1)) {
          if (skip) {
            skip <- FALSE
            next
          }
          w <- which(vec[i] == first_symbols & vec[i+1] == second_symbols)
          if (length(w) == 0) {
            result[res_ind] <- vec[i]
          } else {
            result[res_ind] <- new_symbols[w]
            skip <- TRUE
          }
          res_ind <- res_ind + 1
        }
        if (!skip) {
          result[res_ind] <- vec[n]
          res_ind <- res_ind + 1
        }
        result[seq_len(res_ind-1)]
      })
    },
    .compress = function(matrix, re_pair_method, verbose) {
      if (verbose) message("################## Grammar Compression ##################")
      diff_index_list <- private$.to_diff_index_list(matrix)
      max_symbol <- max(unlist(diff_index_list))
      new_symbol <- max_symbol + 1L
      if (verbose) progress_bar <- txtProgressBar(0, max_symbol, style=3)
      counting <- switch(re_pair_method, lossy=lossy_counting, freq=freq_counting, lossless_counting)
      k <- switch(re_pair_method, plain=1, max_symbol)
      while(TRUE) {
        counts <- counting(diff_index_list, max_hash_size=private$.ncol, l = private$.nrow)
        counts <- rev(sort(counts))
        max_pair <- counts[1]
        if (verbose) setTxtProgressBar(progress_bar, max_symbol - max_pair)
        if (max_pair == 1) break
        k_d <- min(k, sum(counts > 1))
        replace_pairs <- names(counts[seq_len(k_d)])
        diff_index_list <- private$.replace(diff_index_list, replace_pairs, new_symbol:(new_symbol+k_d-1))
        new_symbol <- new_symbol + k_d
      }
      private$.compressed_matrix <- diff_index_list
      private$.grammer_trees <- setNames(names(private$.hash_table), private$.hash_table)
      private$.terminal_symbols <- P <- seq_len(max_symbol)
      names(P) <- P
      for (i in names(private$.grammer_trees)) {
        P[i] <- sum(P[decode_pair(private$.grammer_trees[i])])
      }
      private$.P <- P
      private$.hash_table <- character()
      if (verbose) setTxtProgressBar(progress_bar, max_symbol)
    },
    .recursion = function(i, j, Z, u) {
      if (Z %in% private$.terminal_symbols) {
        if (u + Z == j) {
          private$.R <- c(private$.R, i)
        }
        return(NULL)
      }
      if (u + private$.P[as.character(Z)] == j) {
        private$.R <- c(private$.R, i)
        return(NULL)
      }
      ZLR <- decode_pair(private$.grammer_trees[as.character(Z)])
      left_score <- u + private$.P[as.character(ZLR[1])]
      if (left_score == j) {
        private$.R <- c(private$.R, i)
        return(NULL)
      } else if (left_score > j) {
        private$.recursion(i, j, ZLR[1], u)
      } else {
        private$.recursion(i, j, ZLR[2], left_score)
      }
    }
  ),
  public = list(
    initialize = function(matrix, re_pair_method, verbose = TRUE) {
      private$.ncol <- ncol(matrix)
      private$.nrow <- nrow(matrix)
      private$.compress(matrix, re_pair_method = re_pair_method, verbose = verbose)
    },
    access_row = function(index) {
      diff_indexs <- as.character(private$.compressed_matrix[[index]])
      uncompress <- function(diff_indexs) {
        diff_indexs <- as.character(diff_indexs)
        unlist(lapply(diff_indexs, function(x) {
          if (x %in% private$.terminal_symbols) {
            x
          } else {
            uncompress(decode_pair(private$.grammer_trees[x]))
          }
        }))
      }
      indexs <- cumsum(uncompress(diff_indexs))
      row <- integer(private$.ncol)
      row[indexs] <- 1L
      row
    },
    access_column = function(index) {
      private$.R <- c()
      for (i in seq_len(private$.nrow)) {
        row <- private$.compressed_matrix[[i]]
        u <- c(0L, cumsum(private$.P[as.character(row)]))
        if (any(u == index)) {
          private$.R <- c(private$.R, i)
        } else {
          for (q in seq_along(row)) {
            if (index <= u[q+1]) {
              private$.recursion(i, index, row[q], u[q])
              break
            }
          }
        }
      }
      column <- integer(private$.nrow)
      column[private$.R] <- 1L
      column
    },
    prod = function(matrix) {
      matrix <- as.matrix(matrix)
      result <- matrix(ncol = ncol(matrix), nrow = private$.nrow)
      for (i in seq_len(private$.nrow)) {
        result[i,] <- t(self$access_row(i)) %*% matrix
      }
      result
    },
    tprod = function(matrix) {
      matrix <- as.matrix(matrix)
      result <- matrix(ncol = ncol(matrix), nrow = private$.ncol)
      for (i in seq_len(private$.ncol)) {
        result[i,] <- t(self$access_column(i)) %*% matrix
      }
      result
    }
  )
)

encode_pairs <- function(pairs) {
  paste(pairs[-length(pairs)], pairs[-1], sep="-")
}

decode_pairs <- function(pairs_str) {
  lapply(strsplit(pairs_str, "-"), as.integer)
}

decode_pair <- function(pair_str) {
  decode_pairs(pair_str)[[1]]
}

lossless_counting <- function(diff_index_list, ...) {
  table(Reduce(c, lapply(diff_index_list, encode_pairs)))
}

lossy_counting <- function(diff_index_list, l=length(diff_index_list), ...) {
  pairs <- unlist(lapply(diff_index_list, encode_pairs))

  if (length(pairs) <= length(diff_index_list) * 2) return(lossless_counting(diff_index_list))

  N <- 0
  H <- new.env(hash = TRUE, size = length(pairs), parent = emptyenv())
  delta <- 0
  for (pair_name in pairs) {
    N <- N + 1
    if (exists(pair_name, envir = H)) {
      assign(pair_name, get(pair_name, envir = H) + 1, envir = H)
    } else {
      assign(pair_name, delta + 1, envir = H)
    }
    if (floor(N / l) != delta) {
      delta <- floor(N / l)
      H <- as.list(H)
      H[H < delta] <- NULL
      H <- as.environment(H)
    }
  }
  unlist(as.list(H))
}

freq_counting <- function(diff_index_list, max_hash_size, vacancy_rate = 0.3, ...) {
  pairs <- unlist(lapply(diff_index_list, encode_pairs))
  H <- new.env(hash = TRUE, size = length(pairs), parent = emptyenv())
  for (pair_name in pairs) {
    if (exists(pair_name, envir = H)) {
      assign(pair_name, get(pair_name, envir = H) + 1, envir = H)
    } else {
      hash_size <- length(ls(envir = H))
      if (hash_size >= max_hash_size) {
        while (max_hash_size * (1 - vacancy_rate) < hash_size) {
          H <- lapply(as.list(H), function(x) x-1)
          H[H == 0] <- NULL
          H <- as.environment(H)
          hash_size <- length(ls(envir = H))
        }
      } else {
        assign(pair_name, 1, envir = H)
      }
    }
  }
  unlist(as.list(H))
}

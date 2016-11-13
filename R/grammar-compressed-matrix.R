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
      index_list <- apply(data, 1, function(row) which(row == 1))
      diff_index_list <- lapply(index_list, function(vec) {
        if (length(vec) == 0) {
          integer(0)
        } else {
          c(vec[1], diff(vec))
        }
      })
    },
    .encode_pairs = function(pairs) paste(head(pairs, -1), tail(pairs, -1), sep="-"),
    .decode_pairs = function(pairs_str) lapply(strsplit(pairs_str, "-"), as.integer),
    .decode_pair = function(pair_str) private$.decode_pairs(pair_str)[[1]],
    .lossless_counting = function(diff_index_list) {
      table(Reduce(c, lapply(diff_index_list, private$.encode_pairs)))
    },
    .replace = function(diff_index_list, replace_pair_str, new_symbol) {
      private$.hash_table[replace_pair_str] <- new_symbol
      replaces <- private$.decode_pair(replace_pair_str)
      lapply(diff_index_list, function(vec) {
        n <- length(vec)
        result <- integer(0)
        i <- 1
        while (i <= n) {
          if (i < n && vec[i] == replaces[1] && vec[i+1] == replaces[2]) {
            result <- c(result, new_symbol)
            i <- i + 1
          } else {
            result <- c(result, vec[i])
          }
          i <- i + 1
        }
        Filter(function(x) x >= 0, result)
      })
    },
    .compress = function(matrix) {
      diff_index_list <- private$.to_diff_index_list(matrix)
      max_symbol <- max(Reduce(c, diff_index_list))
      new_symbol <- max_symbol + 1L
      while(TRUE) {
        counts <- private$.lossless_counting(diff_index_list)
        max_pair <- counts[which.max(counts)]
        cat(sprintf("%d.", max_pair))
        if (max_pair == 1) break
        replace_pair <- names(max_pair)
        diff_index_list <- private$.replace(diff_index_list, replace_pair, new_symbol)
        new_symbol <- new_symbol + 1
      }
      private$.compressed_matrix <- diff_index_list
      private$.grammer_trees <- setNames(names(private$.hash_table), private$.hash_table)
      private$.terminal_symbols <- P <- seq_len(max_symbol)
      names(P) <- P
      for (i in names(private$.grammer_trees)) {
        P[i] <- sum(P[private$.decode_pair(private$.grammer_trees[i])])
      }
      private$.P <- P
      private$.hash_table <- character()
    },
    .recursion = function(i, j, Z, u) {
      if (Z %in% private$.terminal_symbols) {
        if (u + Z == j) {
          private$.R <- c(private$.R, i)
        }
        return(NULL)
      }
      ZLR <- private$.decode_pair(private$.grammer_trees[as.character(Z)])
      if (u + private$.P[as.character(ZLR[1])] <= j) {
        private$.recursion(i, j, ZLR[1], u)
      }
      private$.recursion(i, j, ZLR[2], u + private$.P[as.character(ZLR[1])])
    }
  ),
  public = list(
    initialize = function(matrix) {
      private$.ncol <- ncol(matrix)
      private$.nrow <- nrow(matrix)
      private$.compress(matrix)
    },
    access_row = function(index) {
      diff_indexs <- as.character(private$.compressed_matrix[[index]])
      uncompress <- function(diff_indexs) {
        diff_indexs <- as.character(diff_indexs)
        unlist(lapply(diff_indexs, function(x) {
          if (x %in% names(private$.grammer_trees)) {
            uncompress(private$.decode_pair(private$.grammer_trees[x]))
          } else {
            x
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
        for (q in seq_along(row)) {
          if (index <= u[q+1]) {
            private$.recursion(i, index, row[q], u[q])
            break
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

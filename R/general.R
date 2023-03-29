update_renv <- function(){renv::purge("lunafuns");renv::hydrate()}

#' @export ylppa
ylppa <- function(x, ...) {
  dots <- list(...)
  map(dots, \(z)if(!is.function(z))stop("... must only contain functions"))
  setNames(map(dots, \(z)z(x)),
           sapply(as.list(substitute(list(...)))[-1], deparse))
}

#' @export transbind
transbind <- function(...)purrr::map(purrr::list_transpose(list(...)), dplyr::bind_rows)

#' @export parallel_lists
parallel_lists <- function(..., n_args = 2) {
  dots <- substitute(list(...))[-1]
  n <- length(dots)
  if(!identical(n%%n_args, 0))stop("dots length must be a multiple of n_args")
  pos <- seq_len(n/n_args)
  out <- list()
  length(out) <- n_args

  for(i in pos) {
    for(j in seq_len(n_args)) {
      out[[j]][[i]] <- dots[[n_args*(i-1) + j]]
    }
  }

  names(out) <- paste0("arg", seq_len(n_args))
  out
}

#' @export `%p0%`
`%p0%` <- function(x, y){paste0(x, y)}

#' @export ford
ford <- function(...){
  x <- c(...)
  factor(x, levels = unique(x), ordered = TRUE)
}

#' Get a quoted list
#'
#' Designed for use with `targets::tar_map()`, see https://github.com/ropensci/tarchetypes/discussions/105
#'
#' @param ... List contents, can be named
#'
#' @return An expression to create a list with the given contents
#' @export
#'
#' @examples
#'  expr1 <- qist(arg1 = 123, arg2 = letters[1:3], "done")
#'  expr2 <- quote(list(arg1 = 123, arg2 = letters[1:3], "done"))
#'  identical(expr1,expr2)
qist <- function(...){substitute(list(...))}

#' Get indexes matching `y` with each value of `x`
#'
#' @param x,y Some vectors
#'
#' @return A vector of indices of elements in `y` that match `x`, in the order given by `x`
#' @export
#'
#' @examples
#' set.seed(1)
#' which_order(c("b", "d"), sample(letters[1:4], 10, replace = TRUE))
which_order <- function(x, y) {
  as.vector(unlist(apply(sapply(x, `==`, y), 2, which)))
}

#' @export check_length
check_length <- function(x, ...) {
  if(!length(x)%in%c(...))stop("Incompatible length for " %p0%
                                 deparse(substitute(x)))
}

#' Get the names of all LHS variables
#'
#' @param formula_list A list of formulas
#'
#' @return A character vector with the names of unique variables on each LHS, in order of appearance.
#' @export
#'
#' @examples get_lhs_vars(list(y~x, z~y, w+z~x))
get_lhs_vars <- function(formula_list) {
  lhs <- character()
  for(i in formula_list)lhs <- c(lhs, all.vars(i[[2]]))
  unique(lhs)
}

#' Indicate which formulas contain `|mi()` on the LHS
#'
#' @param formula_list A list of formulas
#'
#' @return A logical vector that returns `TRUE` where `|mi()` is found in the formula's LHS.
#' @export
#'
#' @examples get_lhs_mi(list(y ~ x, z | mi() ~ y, w ~ mi(y)))
get_lhs_mi <- function(formula_list) {
  has_mi <- logical()
  for(i in formula_list){
    has_mi <- c(has_mi, all(c("|", "mi()")%in%as.character(i[[2]])))
  }
  has_mi
}

#' Create a dependency matrix from a list of formulas
#'
#' @param formula_list A list containing `formula` class objects
#'
#' @return A square matrix, with names and size equal to the unique variables identified in the formulas. For a given variable's row, values of `1` indicate the variables that depend it. For a given variable's column, values of `1` indicate the variables it depends on.
#' @export
#'
#' @examples
#' {
#'   testf <- list(x11 ~ y1, x12 ~ y1, x13 ~ y1,
#'                 x21 ~ y2, x22 ~ y2, x23 ~ y2,
#'                 x31 ~ y3, x32 ~ y3, x33 ~ y3,
#'                  y3 ~ y1, y2 ~ y1 + z21 + z22)
#'   formlst_depmat(testf)
#' }
formlst_depmat <- function(formula_list) {
  lhs <- get_lhs_vars(formula_list)

  dep_mat<- matrix(0, nrow = length(lhs), ncol = length(lhs),
                   dimnames = list(lhs, lhs))

  for(v in lhs) {
    deps <- rep(0, length(lhs))
    for(i in formula_list) {
      if(all.vars(i[[2]]) == v) {
        dep_i <- which(lhs %in% all.vars(i[[3]]))
        deps[dep_i] <- deps[dep_i] + 1
      }
    }
    dep_mat[,v] <- deps
  }
  dep_mat
}

#' Create a list of variables ordered according to their dependency matrix
#'
#' @param dep_mat A dependency matrix made with `fomrlst_depmat`
#'
#' @return A list where each element contains a vector of variable names that depend on variables in the previous level. The first element holds variables with no dependencies.
#' @export
#'
#' @examples
#' {
#'   testf <- list(x11 ~ y1, x12 ~ y1, x13 ~ y1,
#'                 x21 ~ y2, x22 ~ y2, x23 ~ y2,
#'                 x31 ~ y3, x32 ~ y3, x33 ~ y3,
#'                  y3 ~ y1, y2 ~ y1 + z21 + z22)
#'   depmat_ordlist(formlst_depmat(testf))
#' }
depmat_ordlist <- function(dep_mat) {
  ord_list <- list(names(which(colSums(dep_mat) == 0)))

  while(sum(dep_mat)>0) {
    v_left <- setdiff(colnames(dep_mat), unlist(ord_list))
    dep_mat <- dep_mat[v_left, v_left]

    ord_list <- c(ord_list, list(names(which(colSums(dep_mat) == 0))))
  }
  ord_list
}


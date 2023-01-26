#' @export `%p0%`
`%p0%` <- function(x, y){paste0(x, y)}

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
  lhs_vars <- character()
  for(i in formula_list) {
    lhs_vars <- c(lhs_vars, all.vars(i[[2]]))
  }
  lhs_vars <- unique(lhs_vars)

  dep_mat<- matrix(0, nrow = length(lhs_vars), ncol = length(lhs_vars),
                   dimnames = list(lhs_vars, lhs_vars))

  for(v in lhs_vars) {
    deps <- rep(0, length(lhs_vars))
    for(i in formula_list) {
      if(all.vars(i[[2]]) == v) {
        dep_i <- which(lhs_vars %in% all.vars(i[[3]]))
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

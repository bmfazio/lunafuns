#' Make a function to generate argument tables
#'
#' @param ... Argument names for the desired function, possibly including default values
#'
#' @return A function taking only the arguments passed as `...`, which generates a tibble with the argument names as columns.
#' @export create_arg_maker
#'
#' @examples create_arg_maker(N, M, mu = 0, sigma = 1)
create_arg_maker <- function(...) {
  funcstr <- paste0("function", substring(deparse(substitute(list(...))), 5))
  emptyfun <- paste0(funcstr,"{}")
  tryCatch(eval(parse(text = emptyfun)),
           error = function(x)stop("Invalid expression for function arguments.
                                   Ensure no unnamed values are passed.",
                                   call. = FALSE))
  args <- formals(eval(parse(text = emptyfun)))
  if(any(sapply(args, length)!=1)){
    warning("Not all arguments are length one, may cause unexpected results.")
  }
  out <- eval(parse(text=paste0(funcstr,"ezlist_tibble(tibble::tibble(",
                                paste(names(args),collapse=", "),"))")))
}

#' @export ezlist_tibble
ezlist_tibble <- function(x)dplyr::mutate_all(x, lunafuns_ezlist)

#' @export make_filterable
make_filterable <- function(..., filter_col) {
  if(!is.character(filter_col)|length(filter_col)!=1){
    stop("filter_col must be a length-1 string")
  }
  tb <- dplyr::bind_rows(...)
  if(!filter_col%in%colnames(tb)){
    stop("filter_col must refer to a column in the data frames passed")
  }
  tb[,filter_col] <- sapply(tb[,filter_col], as.character)

  function(filter_string = NULL, exact = FALSE) {
    if(is.null(filter_string))return(tb)
    if(exact)return(dplyr::filter(tb, !!rlang::sym(filter_col) == filter_string))
    dplyr::filter(tb, grepl(filter_string, !!rlang::sym(filter_col)))
  }
}

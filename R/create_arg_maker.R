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
           error = function(x)stop("Invalid expression for function arguments"))
  args <- formals(eval(parse(text = emptyfun)))
  if(any(sapply(args, length)!=1)){
    warning("Not all arguments are length one, may cause unexpected results")
  }
  out <- eval(parse(text=paste0(funcstr,"ezlist_tibble(tibble::tibble(",
                                paste(names(args),collapse=", "),"))")))
}

#' @export ezlist_tibble
ezlist_tibble <- function(x)dplyr::mutate_all(x, lunafuns_ezlist)

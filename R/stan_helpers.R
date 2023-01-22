#' Compile Stan model to desired location
#'
#' @param source Path to a Stan model file
#' @param where Folder where the compiled executable will be placed
#'
#' @return The relative path of the compiled executable
#' @export
sh_compile_path <- function(source, where = ".") {
  model_obj <- cmdstanr::cmdstan_model(source, compile = FALSE)
  model_obj$compile(dir = where)
  file.path(where, fs::path_ext_remove(basename(source)))
}

#' @export targetsh_make_struct
targetsh_make_struct <- function(x)as.list(parse(text = sapply(x, function(y)paste(capture.output(dput(y)), collapse = ""))))

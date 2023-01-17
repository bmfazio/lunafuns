#' @export
new_lunafuns_ezlist <- function(x = list()) {
  x <- purrr::map(x, identity)
  vctrs::vec_assert(x, list())
  vctrs::new_vctr(x, class = "lunafuns_ezlist")
}
#' @export
lunafuns_ezlist <- function(x = list()) {
  new_lunafuns_ezlist(x)
}

#' @importFrom vctrs vec_ptype2 vec_cast
NULL
#' @export
vec_ptype2.lunafuns_ezlist.lunafuns_ezlist <- function(x, y, ...) x
#' @export
vec_ptype2.lunafuns_ezlist.logical <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.logical.lunafuns_ezlist <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.lunafuns_ezlist.integer <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.integer.lunafuns_ezlist <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.lunafuns_ezlist.double <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.double.lunafuns_ezlist <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.lunafuns_ezlist.character <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.character.lunafuns_ezlist <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.lunafuns_ezlist.list <- function(x, y, ...) lunafuns_ezlist()
#' @export
vec_ptype2.list.lunafuns_ezlist <- function(x, y, ...) lunafuns_ezlist()

#' @export
vec_cast.lunafuns_ezlist.lunafuns_ezlist <- function(x, to, ...) x
#' @export
vec_cast.lunafuns_ezlist.logical <- function(x, to, ...) lunafuns_ezlist(x)
#' @export
vec_cast.lunafuns_ezlist.integer <- function(x, to, ...) lunafuns_ezlist(x)
#' @export
vec_cast.lunafuns_ezlist.double <- function(x, to, ...) lunafuns_ezlist(x)
#' @export
vec_cast.lunafuns_ezlist.character <- function(x, to, ...) lunafuns_ezlist(x)
#' @export
vec_cast.lunafuns_ezlist.list <- function(x, to, ...) lunafuns_ezlist(x)

#' @importFrom pillar pillar_shaft
NULL
#' @export
pillar_shaft.lunafuns_ezlist <- function(x, ...) {
  full_desc <- purrr::map(
    x, function(y)paste(capture.output(invisible(dput(y))), collapse=""))

  pillar::new_pillar_shaft_simple(full_desc,width = 30, min_width = 8,
                                   type_sum = "ezlist")
}

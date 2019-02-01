#' Create a transition matrix
#' --------------------------
#'
#' In a competing risk setting state 2 is reachable from state 1.
#' State 3 is absorbing, reachable from states 1 and 2.
#' @param  x:  List of possible transitions.
#'     x[[i]] consists of a vector of state numbers
#'     reachable from state i.
#' @param names: Character vector of state names, having the same length
#'        as x.
#' @examples
#' \donttest{
#'    # For a competing risks setting
#'    trans_mat( x = list( c(2, 3), c(3), c() ) )
#'    # For consecutive states
#'    trans_mat( x = list(c(2), 3, c() ) )
#'    }
#' @references
#' Liesbeth C. de Wreede, Marta Fiocco, Hein Putter
#' (2011). mstate: An R Package for the Analysis of
#' Competing Risks and Multi-State Models. Journal of
#' Statistical Software, 38(7), 1-30. URL
#' http://www.jstatsoft.org/v38/i07/.

trans_mat <- function(x, names) {

  if ( !is.list(x) )
    stop2("x must be a list")

  ns <- length(x) ## number of states
  tmat <- matrix(NA, nrow = ns, ncol = ns) ## transition matrix

  if ( missing(names) ) {
    if ( !is_null( base::names(x) ) ) {
      namesList <- list(from = base::names(x), to = base::names(x))
    } else {
      namesList <- list(from = paste("State", seq(nrow(tmat))),
                        to = paste("State", seq(nrow(tmat))))
    }
  } else {
    if ( length(names) != ns )
      stop2("length of 'names' must equal length of 'x'")
    namesList <- list(from = names, to = names)
  }
  idxmat <- cbind(unlist(lapply(seq(ns),
                                function(i, y){
                                  rep(i, length(y[[i]]))}, x)),
                  unlist(x))
  if ( max(idxmat) > ns )
    stop2("Largest state in transition list exceeds number of states")
  tmat[idxmat] <- seq(nrow(idxmat))
  dimnames(tmat) <- namesList
  tmat
}

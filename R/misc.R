#------------ General
# Create a named list using specified names or, if names are omitted, using the
# names of the objects in the list
#
# @param ... Objects to include in the list.
# @return A named list.
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }

  return(out)
}

# Stop without printing call
stop2    <- function(...) stop(..., call. = FALSE)

# Immediate warning without printing call
warning2 <- function(...) warning(..., immediate. = TRUE, call. = FALSE)

# Shorthand for suppress warnings
SW <- function(expr) base::suppressWarnings(expr)

# Check if an object is NULL
is_null <- function(x) {
  is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)
}

# Safe deparse
safe_deparse <- function(expr) deparse(expr, 500L)

# -------------- Sanity checks ---------#

# Check that a tmat is concordant with
# a formula list
check_trans <- function(f, t){
  if(length(f) != sum(!is.na(t)) ){
    stop2("Formula list has to include a formula for each transition as expressed in the transition matrix")
  }
}


# ------------- Helpers ---------------#
#Compute point estimates and standard errors from pointwise vectors
#
# @param x A matrix.
# @return An ncol(x) by 2 matrix with columns 'Estimate' and 'SE'
#   and rownames equal to colnames(x).
#
table_of_estimates <- function(x) {
  out <- cbind(
    Estimate = colSums2(x),
    SE = sqrt(nrow(x) * colVars(x))
  )
  rownames(out) <- colnames(x)
  return(out)
}


# Check formula object
#
# @param formula The user input to the formula argument.
# @param needs_response A logical; if TRUE then formula must contain a LHS.
validate_formula <- function(formula, needs_response = TRUE) {

  if (!inherits(formula, "formula")) {
    stop2("'formula' must be a formula.")
  }

  if (needs_response) {
    len <- length(formula)
    if (len < 3) {
      stop2("'formula' must contain a response.")
    }
  }
  as.formula(formula)
}

# Check object is a Surv object with a valid type
#
# @param x A Surv object; the LHS of a formula evaluated in a data frame environment.
# @param ok_types A character vector giving the allowed types of Surv object.
validate_surv <- function(x, ok_types = c("right", "counting",
                                          "interval", "interval2")) {
  if (!inherits(x, "Surv"))
    stop2("LHS of 'formula' must be a 'Surv' object.")
  if (!attr(x, "type") %in% ok_types)
    stop2("Surv object type must be one of: ", comma(ok_types))
  x
}



# ------- Helpers data ------------------ #
# Parse the model formula
#
# @param formula The user input to the formula argument.
# @param data The user input to the data argument (i.e. a data frame).
parse_formula <- function(formula, data) {

  formula <- validate_formula(formula, needs_response = TRUE)

  lhs      <- lhs(formula) # full LHS of formula
  lhs_form <- reformulate_lhs(lhs)

  rhs        <- rhs(formula)         # RHS as expression
  rhs_form   <- reformulate_rhs(rhs) # RHS as formula
  rhs_terms  <- terms(rhs_form, specials = "tde")
  rhs_vars   <- rownames(attr(rhs_terms, "factors"))

  allvars <- all.vars(formula)
  allvars_form <- reformulate(allvars)

  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")

  if (type == "right") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
    min_t    <- 0
    max_t    <- max(surv[, "time"])
  } else if (type == "counting") {
    tvar_beg <- as.character(lhs[[2L]])
    tvar_end <- as.character(lhs[[3L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- min(surv[, "start"])
    max_t    <- max(surv[, "stop"])
  } else if (type == "interval") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- 0
    max_t    <-  max(surv[, c("time1", "time2")])
  } else if (type == "interval2") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
    min_t    <- 0
    max_t    <- max(surv[, c("time1", "time2")])
  }

  sel <- attr(rhs_terms, "specials")$tde

  if (!is.null(sel)) { # model has tde

    # replace 'tde(x, ...)' in formula with 'x'
    tde_oldvars <- rhs_vars
    tde_newvars <- sapply(tde_oldvars, function(oldvar) {
      if (oldvar %in% rhs_vars[sel]) {
        tde <- function(newvar, ...) { # define tde function locally
          safe_deparse(substitute(newvar))
        }
        eval(parse(text = oldvar))
      } else oldvar
    }, USE.NAMES = FALSE)
    term_labels <- attr(rhs_terms, "term.labels")
    for (i in sel) {
      sel_terms <- which(attr(rhs_terms, "factors")[i, ] > 0)
      for (j in sel_terms) {
        term_labels[j] <- gsub(tde_oldvars[i],
                               tde_newvars[i],
                               term_labels[j],
                               fixed = TRUE)
      }
    }
    tf_form <- reformulate(term_labels, response = lhs)

    # extract 'tde(x, ...)' from formula and construct 'bs(times, ...)'
    tde_terms <- lapply(rhs_vars[sel], function(x) {
      tde <- function(vn, ...) { # define tde function locally
        dots <- list(...)
        ok_args <- c("df")
        if (!isTRUE(all(names(dots) %in% ok_args)))
          stop2("Invalid argument to 'tde' function. ",
                "Valid arguments are: ", comma(ok_args))
        df <- if (is.null(dots$df)) 3 else dots$df
        degree <- 3
        if (df == 3) {
          dots[["knots"]] <- numeric(0)
        } else {
          dx <- (max_t - min_t) / (df - degree + 1)
          dots[["knots"]] <- seq(min_t + dx, max_t - dx, dx)
        }
        dots[["Boundary.knots"]] <- c(min_t, max_t)
        sub("^list\\(", "bs\\(times__, ", safe_deparse(dots))
      }
      tde_calls <- eval(parse(text = x))
      sel_terms <- which(attr(rhs_terms, "factors")[x, ] > 0)
      new_calls <- sapply(seq_along(sel_terms), function(j) {
        paste0(term_labels[sel_terms[j]], ":", tde_calls)
      })
      nlist(tde_calls, new_calls)
    })
    td_basis <- fetch(tde_terms, "tde_calls")
    new_calls <- fetch_(tde_terms, "new_calls")
    td_form <- reformulate(new_calls, response = NULL, intercept = FALSE)

  } else { # model doesn't have tde
    tf_form  <- formula
    td_form  <- NULL
    td_basis <- NULL
  }

  nlist(formula,
        lhs,
        rhs,
        lhs_form,
        rhs_form,
        tf_form,
        td_form,
        td_basis,
        fe_form = rhs_form, # no re terms accommodated yet
        re_form = NULL,     # no re terms accommodated yet
        allvars,
        allvars_form,
        tvar_beg,
        tvar_end,
        dvar,
        surv_type = attr(surv, "type"))
}


# Extract LHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
lhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[2L]]
  } else {
    out <- NULL
  }
  out
}


# Find starting
#
# @param a transition matrix
# @param a transition
match_starting <- function(t, m){
  which(apply(m, 1, function(x) return(t %in% x )))
}

# Find next state
#
# @param a transition matrix
# @param a transition
match_to <- function(t, m){
  which(apply(m, 2, function(x) return(t %in% x )) )

}

# Find competing state
#
# @param a transition matrix
# @param a initial state
# @param a ending state
match_competing <- function(m, r, t){
  m <- m[r, ][-match(names(t), names(m[r,]))]
  m[!is.na(m)]
}

# Reformulate as LHS of a formula
#
# @param x A character string or expression object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_lhs <- function(x) {
  #x <- deparse(x, 500L)
  x <- formula(substitute(LHS ~ 1, list(LHS = x)))
  x
}

# Extract RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
rhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[3L]]
  } else {
    out <- x[[2L]]
  }
  out
}


# Return a data frame with NAs excluded
#
# @param formula The parsed model formula.
# @param data The user specified data frame.
make_model_data2 <- function(data, formula, cens_formula ) {

  time_var <- data[ ,formula$tvar_end]
  status_var <- data[ ,formula$dvar]

  if(length(cens_formula)){
    censor_time <- lapply(cens_formula, function(c)
      data[ ,c$tvar_end] )
    for(i in seq_along(censor_time)){

    }

  }

  data <- data[data[aux_formula$dvar] == aux_cens, ]
  data <- data[data[formula$tvar_end] > 0, ] # remove 0 time
  data <- data[data[formula$dvar] == cens, ]
  mf <- model.frame(formula$tf_form, data, na.action = na.pass)
  include <- apply(mf, 1L, function(row) !any(is.na(row)))
  data[include, , drop = FALSE]
}




# # Return a data frame with NAs excluded
# #
# # @param formula The parsed model formula.
# # @param data The user specified data frame.
# make_model_data <- function(formula, aux_formula, data, cens, aux_cens ) {
#
#   data <- data[data[aux_formula$dvar] == aux_cens, ]
#   data <- data[data[formula$tvar_end] > 0, ] # remove 0 time
#   data <- data[data[formula$dvar] == cens, ]
#   mf <- model.frame(formula$tf_form, data, na.action = na.pass)
#   include <- apply(mf, 1L, function(row) !any(is.na(row)))
#
#
#   data[include, , drop = FALSE]
# }

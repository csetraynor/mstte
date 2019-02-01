ms_prep <- function(formula,
                    data,
                    tm,
                    recovery,
                    stateInit,
                    ids,
                    status,
                    times,
                    keep){

  if(missing(data)){
    stop2("A dataset has to be provided.")
  }
  if(missing(tm)){
    stop2("A transition matrix has to be provided. One can be created e.g. via trans_mat()")
  }
  if(missing(stateInit)){
    stateInit <- NULL
  }
  ns <- max(tm, na.rm = TRUE)
  if(missing(formula)){
    formula <- NULL
  }
  if(missing(recovery)){
    recovery = FALSE
  }
  if(missing(ids)){
    ids <- seq_len(nrow(data))
  }

  check_trans(tm, formula)


  states.to <- apply(!is.na(tm),1,sum)
  absorbing <- which(states.to==0)
  states.from <- apply(!is.na(tm),2,sum)
  if(!recovery){
    startings <- which(states.from==0)
  } else if(!is_null(stateInit)){
    if(is.character(stateInit)){
      stateInit <- match(stateInit, colnames(tm))
    }
    startings <- which(is.na(tm[1,])[stateInit])
  } else {
    warning2("Assuming a first initial state from the first column of the transition matrix")
    startings <- which(is.na(tm[1,])[1])
  }


  if(!is.null(formula)){

    # check_trans
    check_trans(tm, formula)

    # parse formula
    formula <- lapply(formula, function(f) parse_formula(f, data))


    data <- lapply(seq_len(ns), function(s){
      nr <- match_starting(s, tm)
      nt <- match_to(s, tm)
      ni <- match_censors(tm, nr, nt)

      if(recovery){

      } else {
        formula_cens =  formula[ni]
      }

      data <- make_model_data2(data,
                               formula = formula[[s]],
                               formula_cens = )




    })





  }

}



  # Return a data frame with NAs excluded
  #
  # @param formula The parsed model formula.
  # @param data The user specified data frame.
  make_model_data <- function(formula, aux_formula, data, cens, aux_cens ) {

    data <- data[data[aux_formula$dvar] == aux_cens, ]
    data <- data[data[formula$tvar_end] > 0, ] # remove 0 time
    data <- data[data[formula$dvar] == cens, ]
    mf <- model.frame(formula$tf_form, data, na.action = na.pass)
    include <- apply(mf, 1L, function(row) !any(is.na(row)))


    data[include, , drop = FALSE]
  }

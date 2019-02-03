





ms_stan <- function(formula,
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
    ns <- max(tm, na.rm = TRUE)

    # parse formula
    formula <- lapply(formula, function(f) parse_formula(f, data))

    # prepare data
    data <- lapply(seq_len(ns), function(s){
      make_model_data2(data, formula = formula[[s]], s)
    })






  }

}

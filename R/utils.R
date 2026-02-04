#' @keywords internal
#'
scaling_zscore = function(x){
  if (typeof(x) == "double"){

    if(length(x) == 1){ ## in case the vector has length 1 -- return 0
      return(0)
    }
    if(sd(x, na.rm = TRUE) > 0){
      return((x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
    } else{
      return((x - mean(x, na.rm = TRUE)))
    }
  } else {return(x)}
}
#' @keywords internal
#'
scaling_modified_zscore = function (x) {
  if (typeof(x) == "double") {
    if (mad(x, na.rm = TRUE) != 0) {
      return(0.6745 * (x - median(x))/mad(x))
    }
    else {
      return(0.6745 * (x - median(x)))
    }
  }
  else {
    return(x)
  }
}

mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])

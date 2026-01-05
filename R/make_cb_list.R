#' IN-PROGRESS functions that creates a list of cross-basis objects to try
#'
#' @param x the data.frame with the var_list data in it
#' @param var_list
#' @param argvar_list
#' @param lag_list
#' @param arglag_list
#'
#' @returns
#'
#' @examples
make_cb_list <- function(x, var_list, argvar_list, lag_list, arglag_list) {
  #' x:
  #' var_list:    character vector of columns of X that you use, e.g., c('tmean')
  #' argvar_list: list(set1 = list(fun="bs", degree = 2, knots = list(vk)))
  #' lag_list:    c(5, 30)
  #' arglag_list: list(set1 = list(knots = list(lk))), same fcn as argvar_list

  library(dlnm); library(splines); library(tidyverse)

  # --- function to summary arguments ---
  get_summary_string <- function(x) {
    paste(sapply(names(x), function(n) {
      value <- x[[n]]
      ifelse(length(value) > 1,
             paste0(n, "=len[", length(value), "]"),
             paste0(n, "=", ifelse(is.numeric(value), round(value, 3), value)))
    }), collapse = ",")
  }

  # --- x ---
  stopifnot(is.data.frame(x))

  # --- var ---
  stopifnot(typeof(var_list) == 'character')
  stopifnot(all(var_list %in% colnames(x)))
  stopifnot(all(sapply(var_list, function(js) is.numeric(unlist(x[,js])))))
  var_index <- 1:length(var_list)
  var_name <- var_list
  var_list <- lapply(var_list, function(s) unlist(x[, s]))

  # --- lag ---
  stopifnot(is.numeric(lag_list))
  stopifnot(all(lag_list %% 1 == 0))

  # --- argvar ---
  argvar_list_expanded <- expand_arglist(argvar_list)

  # --- arglag ---
  arglag_list_expanded <- expand_arglist(arglag_list)

  # --- combine ---
  cb_df <- expand_grid(var_index, argvar_list_expanded,
                       lag_list, arglag_list_expanded)
  cb_df$var_list <- lapply(var_index, function(i) var_list[[i]])
  cb_df$var_name <- sapply(var_index, function(i) var_name[i])

  # --- pass into crossbasis
  cb_out <- vector('list', nrow(cb_df))
  for(i in 1:nrow(cb_df)) {
    this_x      = cb_df$var_list[i][[1]]
    this_argvar = cb_df$argvar_list_expanded[i][[1]]
    this_lag    = cb_df$lag_list[i][[1]]
    this_arglag = cb_df$arglag_list_expanded[i][[1]]
    # crossbasis
    cb_out[[i]] <- crossbasis(
      x      = this_x,
      argvar = this_argvar,
      lag    = this_lag,
      arglag = this_arglag)
    # add name attribute
    attr(cb_out[[i]], 'var_name') <- cb_df$var_name[i]
    # need something else:
    attr(cb_out[[i]], 'summary_str') <- paste0(
      sprintf("var=%s:{",cb_df$var_name[i]),
      get_summary_string(this_argvar),"}, ",
      sprintf("lag=%i:{", this_lag),
      get_summary_string(this_arglag),"}"
    )
  }

  return(cb_out)

}

#' IN-PROGRESS functions that creates a list of cross-basis objects to try
#'
#' @param arglist
#'
#' @returns
#'
#' @examples
expand_arglist <- function(arglist) {

  # Empty list to hold the output.
  out_l <- list()

  # --- Helper functions ---
  # Define a function to create a grid of all combinations of two inputs.
  expand_grid2 <- function(a, b, b_name) {
    if (nrow(a) == 0) {
      tmp <- NA
      o <- data.frame(tmp)
      o[b_name] <- list(b)
    } else {
      o <- expand_grid(a, b)
      colnames(o)[ncol(o)] <- b_name
    }
    return(o)
  }

  # Define a function to bind two inputs into columns of a data frame.
  cbind2 <- function(a, b, b_name) {
    if (nrow(a) == 0) {
      tmp <- NA
      o <- data.frame(tmp)
    }
    o[b_name] <- list(b)
    return(o)
  }

  # Define a function to unlist the 'knots' elements of a list.
  unlist_knots <- function(lst) {
    for (i in 1:length(lst)) {
      for (name_i in names(lst[[i]])) {
        if (name_i == "knots") {
          stopifnot(length(lst[[i]]["knots"]) == 1)
          lst[[i]]["knots"] <- lst[[i]]["knots"][[1]]
        }
      }
    }
    return(lst)
  }

  # --- Main loop ---
  # Loop over each set of arguments in 'arglist'.
  for (arg_set_i in arglist) {

    # Get the names of the arguments in the current set.
    arg_i_names <- names(arg_set_i)

    # Create an empty data frame to hold the output for the current set.
    out_i <- data.frame()

    # Loop over each argument in the current set.
    for (arg_i in arg_i_names) {

      # Get the value of the current argument.
      y <- arg_set_i[[arg_i]]

      if (arg_i == "knots") {

        # If the argument is a list, expand it into a grid.
        # If it is not a list, bind it into a column of a data frame.
        if (typeof(y) == "list") {
          out_i <- expand_grid2(out_i, y, arg_i)
        } else {
          out_i <- cbind2(out_i, y, arg_i)
        }

      } else {

        # If the argument is not known, issue a warning.
        known_vars <- c("df", "intercept", "degree", "fun")
        if (!(arg_i %in% known_vars)) {
          warning(sprintf(
            'Argument "%s" is not in [%s], behavior unknown.',
            arg_i, paste(known_vars, collapse = ", ")
          ))
        }

        # Assume the argument is a vector, and expand it into a grid.
        y <- arg_set_i[[arg_i]]
        stopifnot(typeof(y) %in% c("double", "character"))
        out_i <- expand_grid2(out_i, y, arg_i)

      }
    }

    # Convert each row of the output data frame into a named list.
    out_l_i <- lapply(1:nrow(out_i), function(i) sapply(out_i[i, ], list))

    # Unlist the 'knots' elements of the output list.
    out_l_i_clean <- unlist_knots(out_l_i)

    # Append the output list for the current set to the overall output list.
    out_l <- append(out_l, out_l_i_clean)
  }

  # return and remove any named $tmp
  lapply(out_l, function(ll) ll[! (names(ll) %in% "tmp")])
}

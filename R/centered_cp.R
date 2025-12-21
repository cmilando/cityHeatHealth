#' Internal function to get centered cp objects
#'
#' Needed a function for this because we do it twice: once for regional RRs
#' and once for BLUP
#'
#' @param argvar
#' @param xcoef
#' @param xvcov
#' @param this_exp
#' @param x_b
#' @param global_cen
#' @param cen
#'
#' @returns
#'
#' @examples
get_centered_cp <- function(argvar, xcoef, xvcov,
                            this_exp, x_b,
                            global_cen, cen) {

  # define grid
  grid <- seq(from =  x_b[1], to = x_b[2], by = 0.1)

  # (1) get onebasis
  basis_x <- do.call("onebasis",
                     modifyList(argvar,
                                list(x = this_exp,
                                     Boundary.knots = x_b)))

  # *********
  # (2) Center basis
  # either MMT or GLOBAL CEN
  # you need boundary knots because the centerpoint is almost
  # always outside of the percentiles in this work
  # so this creates a full range to test over
  if(!is.null(global_cen)) {

    cen = global_cen
    stopifnot(global_cen >= x_b[1] & global_cen <= x_b[2])
    basis_mmt <- do.call("onebasis", modifyList(argvar,
                                                list(x=global_cen,
                                                     Boundary.knots = x_b)))
  } else {

    b2 <- crosspred(basis_x,
                    cen = mean(this_exp),
                    coef = xcoef,
                    vcov = xvcov,
                    model.link = "log",
                    at = grid)

    cen = b2$predvar[which.min(b2$allRRfit)]

    basis_mmt <- do.call("onebasis",
                         modifyList(argvar,
                                    list(x=cen, Boundary.knots = x_b)))
  }

  # *********

  # (3) Center and scale
  basis_cen <- scale(basis_x, center = basis_mmt, scale = FALSE)


  # get the cross-pred object
  # cen is passed forward from before
  # the main reason  you need this for the RR plot
  # and this gives back out BLUP coef and vcov which you can use
  # in the AN calc
  # does it also give back basis_cen?
  centered_cp <- crosspred(basis_cen,
                       cen = cen,
                       coef = xcoef,
                       vcov = xvcov,
                       model.link = "log",
                       at = grid)

  return(list(cp = centered_cp, basis_cen = basis_cen))

}

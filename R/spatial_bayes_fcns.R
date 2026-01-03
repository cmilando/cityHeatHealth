#' Title
#'
#' @param x
#' @param c_id
#' @param comp_ids
#' @param offsets
#'
#' @returns
#' @export
#'
#' @examples
nb2subgraph = function(x, c_id, comp_ids, offsets) {
  # nb2subgraph
  # for a given subcomponent, return graph as lists of node1, node2 pairs
  #
  # inputs:
  # x: nb object
  # c_id: subcomponent id
  # comp_ids: vector of subcomponent ids
  # offsets: vector of subcomponent node numberings
  # returns: list of node1, node2 ids
  #
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        n_links = n_links + length(x[i]);
      }
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        for (j in 1:length(x[[i]])) {
          n2 = unlist(x[[i]][j]);
          if (i < n2) {
            idx = idx + 1;
            node1[idx] = offsets[i];
            node2[idx] = offsets[n2];
          }
        }
      }
    }
  }
  return (list("node1"=node1,"node2"=node2));
}


#' Title
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
nb2graph = function(x) {
  # nb2graph
  #
  # input: nb_object
  # returns: dataframe containing num nodes, num edges,
  #          and a list of graph edges from node1 to node2.
  #
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (x[[i]][1] != 0) {
      n_links = n_links + length(x[[i]]);
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (x[[i]][1] > 0) {
      for (j in 1:length(x[[i]])) {
        n2 = unlist(x[[i]][j]);
        if (i < n2) {
          idx = idx + 1;
          node1[idx] = i;
          node2[idx] = n2;
        }
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}


#' Title
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
indexByComponent = function(x) {
  # indexByComponent
  #
  # input: vector of component ids
  # returns: vector of per-component consecutive node ids
  #
  y = x;
  comps = as.matrix(table(x));
  num_comps = nrow(comps);
  for (i in 1:nrow(comps)) {
    idx = 1;
    rel_idx = 1;
    while (idx <= length(x)) {
      if (x[idx] == i) {
        y[idx] = rel_idx;
        rel_idx = rel_idx + 1;
      }
      idx = idx + 1;
    }
  }
  return(y);
}


#' Title
#'
#' https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/nb_data_funs.R
#'
#' @param x
#'
#' @returns
#' @export
#' @import Matrix
#' @importFrom INLA inla.qinv
#' @examples
scale_nb_components = function(x) {
  N = length(x);
  comp_ids = n.comp.nb(x)[[2]];
  offsets = indexByComponent(comp_ids);

  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  scales = vector("numeric", length=num_comps);
  for (i in 1:num_comps) {
    N_subregions = comps[i,1];
    scales[i] = 0.0;
    if (N_subregions > 1) {
      # get adj matrix for this component
      drops = comp_ids != i;
      nb_tmp = droplinks(x, drops);
      nb_graph = nb2subgraph(nb_tmp, i, comp_ids, offsets);
      adj.matrix = sparseMatrix( i=nb_graph$node1, j=nb_graph$node2, x=1, dims=c(N_subregions,N_subregions), symmetric=TRUE);
      # compute ICAR precision matrix
      Q =  Diagonal(N_subregions, rowSums(adj.matrix)) - adj.matrix;
      # Add a small jitter to the diagonal for numerical stability (optional but recommended)
      Q_pert = Q + Diagonal(N_subregions) * max(diag(Q)) * sqrt(.Machine$double.eps)
      # Compute the diagonal elements of the covariance matrix subject to the
      # constraint that the entries of the ICAR sum to zero.
      Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N_subregions),e=0))
      # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
      scaling_factor = exp(mean(log(diag(Q_inv))))
      scales[i] = scaling_factor;
    }
  }
  return(scales);
}

#' Title
#'
#' @param shp
#' @param ni
#' @param include_self
#'
#' @returns
#' @export
#'
#' @examples
getSW <- function(shp, ni, include_self = T) {


  #
  warning("has to be one polygon per row in `shp`")
  J <- nrow(shp)

  if(ni == 0) {
    return(diag(J))
  }

  #
  nb <- poly2nb(local_shp)
  lb <- nb2listw(nb, style = 'B', zero.policy = TRUE)
  xx <- lb$neighbours

  # get all neighbor lists
  nxi <- nx(xx, ni, include_self)

  # get ones that aren't duplicated
  ## nx <- nxi[which(!duplicated(nxi))]
  # if you don't do this, you can use SW as the Jmatrix
  nx <- nxi

  # turn it into a matrix
  s_rows <- length(nx)

  SW <- matrix(0, nrow  = s_rows, ncol = J)
  for(i in 1:nrow(SW)) {
    SW[i, nx[[i]]] <- 1
  }

  # class(SW) = 'SW'

  return(SW)

}

#' Title
#'
#' @param xx
#' @param ni
#' @param include_self
#'
#' @returns
#' @export
#'
#' @examples
nx <- function(xx, ni, include_self = T) {

  ##
  if(ni == 1) {
    if(include_self) {
      return(lapply(1:length(xx), \(i) sort(c(i, xx[[i]]))))
    } else{
      return(lapply(1:length(xx), \(i) sort(c(xx[[i]]))))
    }

  } else {

    n_out <- vector("list", length(xx))

    for(i in 1:length(xx)) {

      n_minus_1 <- nx(xx, ni = ni - 1)

      s1 <- do.call(c, lapply(n_minus_1[[i]], \(x) n_minus_1[[x]]))

      if(include_self) {
        n_out[[i]] <- sort(unique(c(i, s1)))
      } else {
        n_out[[i]] <- sort(unique(c(s1)))
      }

    }

    return(n_out)

  }

  ##

}


#' Title
#'
#' @param strata_vector
#'
#' @returns
#' @export
#'
#' @examples
getSmat <- function(strata_vector) {
  # create S matrix
  # strata_vector <- data$stratum
  strata_matrix <- matrix(as.integer(strata_vector),
                          nrow = length(strata_vector),
                          ncol = length(strata_vector),
                          byrow = T)

  for(i in 1:length(strata_vector)) {
    strata_matrix[i, ] = 1*(strata_matrix[i, ] == as.integer(strata_vector[i]))
  }

  return(strata_matrix)
}

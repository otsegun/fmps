library(fda)
library(ggplot2)
library(viridis)

## matern function
my.matern <- function(h, b, sigma, nu){
  h[h == 0] <- 1e-10
  num1 <- (sigma ^ 2) * (2 ^ (1 - nu)) / gamma(nu)
  num2 <- (h / b) ^ nu
  num3 <- besselK(x = (h / b), nu = nu)
  return(num1 * num2 * num3)
}

# guassian function
my.Gaussian <- function(h, b, sigma){
  h[h == 0] <- 1e-10
  num1 <- sigma ^ 2 * exp(-(h / b) ^ 2)
  return(num1)
}

# euclidean distance
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))




# simulate matern random field given parameters amd X and Y range
sim_matern_field <- function(x = seq(0, 1, length.out = 70),
                             y = seq(0, 1, length.out = 70),
                             sigma = 1,
                             b = 0.03339183,
                             nu = 0.5){
  ggrid <- expand.grid(x,y)
  stat_cov <- my.matern(h = fields::rdist(ggrid), b = b, nu = nu, sigma = sigma)
  err <- mvtnorm::rmvnorm(n = 1, sigma = stat_cov, method = "chol")
}

# simulate functional random field
# xseq: a vector containing the x range value
# yseq: a vector containing the y range values
# basistype: basis to use for generating functions at each location. Default is B-Spline ("spline")
# tvec: range of t values (the function domain)
# field_params: a list of three vectors which are parameters for generating each
#               random field used in linear combination. Must contain the vectors
#               sigma, beta, nu

simulate_functional_rf <- function(xseq,
                       yseq,
                       basistype = c("spline", "fourier"),
                       nobasis = 10,
                       tvec = seq(from = 0, to = 1, len = 200),
                       field_params) {
  # this model generates spatial functional data with 
  # mean function zero + linear combo of n_basis basis functions
  # and random field. The random field model and its parameters are fixed
  # so that n_basis realizations of this rf models are linearly combined
  # with the basis functions. 
  basistype <- match.arg(basistype)
  if (basistype == "spline") {
    spbasis <- create.bspline.basis(rangeval = c(0, 1), nbasis = nobasis)
  } else if (basistype == "fourier") {
    spbasis <- create.fourier.basis(rangeval = c(0, 1), nbasis = nobasis)
  } else{
    stop("basis currently not supported")
  }
  spbasismatrix <- eval.basis(tvec, spbasis)
  
  # simulate coefficient of basis functions
  sim_data <- list()
  xy <- unname(as.matrix(expand.grid(xseq,yseq)))
  for (i in 1:nobasis) {
    cat("::Generating random field for basis: ", i, "::\n")
    sim_data[[i]] <- sim_matern_field(x = xseq, y = yseq,
                                      sigma = field_params$sigma[i],
                                      b = field_params$beta[i],
                                      nu = field_params$nu[i]) 
    
  }
  
  final_res <- matrix(0, nrow = nrow(xy), ncol = length(tvec))
  
  # linear combination
  for(k in seq(nobasis)) {
    final_res <- final_res + (matrix(sim_data[[k]], ncol =  1) %*% matrix(spbasismatrix[, k], nrow = 1))
  }
  return(list(sim_data = final_res, coords = xy, rndfields = sim_data))
}

#### distances for comparing functional data events

# distance 1
# cand_fdata: candidate functional data event (from training image)
# pfdata: functional data event (from simulation grid)
# wts: weight vector
data_event_dist1 <- function(cand_fdata,
                             pfdata, 
                             wts = NULL) {
  d1 <- dim(cand_fdata)
  d2 <- dim(pfdata)
  difffdata <- cand_fdata - pfdata
  intm_dist <- rowSums((difffdata)^2)
  if(is.null(wts)){
    return(sqrt(mean(intm_dist))) # equal weighting
  }else{
    wts <- wts/(sum(wts))
    return(sqrt(weighted.mean(intm_dist, wts)))
  }
}

# distance 2
# cand_fdata: candidate functional data event (from training image)
# pfdata: functional data event (from simulation grid)
# dmax: maximum distance
data_event_dist2 <- function(cand_fdata,
                             pfdata, 
                             dmax){
  intm_dist <- sqrt(rowSums((cand_fdata - pfdata)^2))
  return(mean(intm_dist/dmax))
  
}
# lag-1 differencing
difference_curves <- function(dt){
  p <- dim(dt)[2]
  dt[,2:p] - dt[, 1:(p-1)]
}

# distance 3
# cand_fdata: candidate functional data event (from training image)
# pfdata: functional data event (from simulation grid)
# wts: weight vector
# deriv_order: order of derivative
data_event_dist_deriv <- function(cand_fdata,
                                  pfdata, 
                                  deriv_order, 
                                  wts = NULL){
  
  for (i in seq(deriv_order)) {
    cand_fdata <- difference_curves(cand_fdata)
    pfdata <- difference_curves(pfdata)
  }
  
  intm_dist <- rowSums((cand_fdata - pfdata)^2)
  
  if(is.null(wts)){
    return(sqrt(mean(intm_dist))) # equal weighting
  }else{
    wts <- wts/(sum(wts))
    return(sqrt(weighted.mean(intm_dist, wts)))
  }
  
}

# distance 4
# cand_fdata: candidate functional data event (from training image)
# pfdata: functional data event (from simulation grid)
# wts: weight vector

data_event_sim_map_muod <- function(cand_fdata, pfdata, wts = NULL){
  dm <- dim(cand_fdata)
  nr <- dm[1]; nc <- dm[2]
  sd_cand <- apply(cand_fdata, 1, sd)
  sd_pfdata <- apply(pfdata, 1, sd)
  mn_cand <- rowMeans(cand_fdata)
  mn_pfdata <- rowMeans(pfdata)
  
  sy <- ay <- my <-  vector(length = nr)
  
  for (qx in 1:nr) {
    cv <- cov(cand_fdata[qx,], pfdata[qx, ])
    sy[qx] <- cv/(sd_cand[qx]*sd_pfdata[qx])
    ay[qx] <- cv/(sd_pfdata[qx]^2)
    my[qx] <- mn_cand[qx] - ay[qx]*mn_pfdata[qx]
  }
  sy <- 1-sy;   ay <- abs(ay-1); my <- abs(my)
  if(is.null(wts)){
    Sy <- mean(sy) 
    Ay <- mean(ay)
    My <- mean(my)
  }else{
    wts <- wts/(sum(wts))
    Sy <- weighted.mean(sy, wts)
    Ay <- weighted.mean(ay, wts)
    My <- weighted.mean(my, wts)
  }
  return(c(Ay, My, Sy))
  
}


# function to convert lag vectors and candidate point index to neighbourhood indices
# used in functional mps
find_cand_ngbs_inds <- function(simgridd, candptind, lagvectors, ccdiff, lenout){
  # convertion of lag vectors and candidate point ind to negbhd inds
  # cand_ngbhd_ind = cand_pt_ind + (len_out * (lag_vector_x2/common_diff)) + (lag_vector_x1/common_diff)
  #                = cand_pt_ind + (lag_vector_x2 * (len_out/common_diff) + lag_vector_x1/common_diff)
  # check to see that simgrid[candidate_ngbs_ind, ] == candidate_ngbs
  nr <- dim(lagvectors)[1]
  candpoint <- simgridd[candptind, ]
  candidate_ngbs <- lagvectors + rep(candpoint, each = nr)
  cand_ngbs_ind <- round(candptind + rowSums(lagvectors * rep(c(1/ccdiff, (lenout/ccdiff)), each = nr)))
  # test
  if(isTRUE(all.equal(simgridd[cand_ngbs_ind, ], candidate_ngbs))){
    return(cand_ngbs_ind)
  }else{
    cat("::Error converting cand pt ind and lag vecs to cand ngbs ind::\n")
    stop("all.equal(simgridd[candidate_ngbs_ind, ], candidate_ngbs) test failed!")
  }
  
}


##### functional mps with gap filling #####
# simgrid: coords of the training image and simulation grid (assumed to be the same)
# tdata: training image with data
# simdata: simulation_grid with some points missing to fill
# missing_points : indices of points to fill in simdata
# nnghbs: number of neighbours
# xseq: range vector of x and y locations (square grids considered)
# dist_metric: the distance metric
# weighted: whether to weight disance metric
# wts_muod: weights for distance 4
fun_mps2 <- function(simgrid,
                    tdata,
                    simdata,
                    missing_points,
                    nnghbs, 
                    xseq, 
                    dist_metric, 
                    weighted = T,
                    wts_muod = c(1/3, 1/3, 1/3) # amp, mag, shape,
                    ){
  cdiff <- xseq[2]-xseq[1]
  len_out <- length(xseq)
  
  dm <- dim(tdata);
  nr <- dm[1];
  nc <- dm[2]
  
  # define random path across simulation grid
  rowss <- 1:nr
  cond_points <- rowss[!(rowss %in%  missing_points)]
  
  if(length(missing_points) > 1){
    simpoints <- sample(missing_points) # shuffle simpoint
  }
  
  lsmp <- length(simpoints)
  
  #find dmax if distance is d2
  if(dist_metric == "d2"){
    tdata_norm <- sqrt(rowSums(tdata^2))
    tnorm_max <- which.max(tdata_norm)
    tnorm_min <- which.min(tdata_norm)
    dmax <- sqrt(sum((tdata[tnorm_max,] - tdata[tnorm_min,])^2))
  }
  
  # start simulation
  counter <- 1
  
  for (i in simpoints) {
    cat("::---Candidate ", counter, "out of ", lsmp, "for location", i, "-::", "\n")
    coords_i <- simgrid[i, ]
    lcond <- length(cond_points)
    
    # find n nearest neighbours to point i
    coords_cond_points <- simgrid[cond_points, ]
    dist_i_cond_points <- apply(coords_cond_points, 1, euc.dist, x2 = coords_i)
    sorted_dists <- sort(dist_i_cond_points, index.return = T)
    ind_nb_pts <- sorted_dists$ix[1:nnghbs] #
    nb_pts <-   cond_points[ind_nb_pts]
    nb_pts_coords <- simgrid[nb_pts, ]
    
    # define search window in TI
    lag_vectors <- nb_pts_coords - rep(coords_i, each = nnghbs)
    if(isTRUE(weighted)) lag_vector_norms <- sqrt(rowSums(lag_vectors^2))
    xy_max_limit <- round(max(xseq) - matrixStats::colMaxs(lag_vectors), 2)
    xy_min_limit <- round(min(xseq) - matrixStats::colMins(lag_vectors),2)
    test_cond <- (simgrid[, 1] <= xy_max_limit[1]) & # problem of rounding
      (simgrid[, 1] >= xy_min_limit[1]) &
      (simgrid[, 2] <= xy_max_limit[2]) &
      (simgrid[, 2] >= xy_min_limit[2])
    search_window_ind <- which(test_cond) # indexes simgrid
    lswi <- length(search_window_ind)
    if(lswi > 1){ # sampling behaves differently when it receives an integer
      search_window_ind <- sample(search_window_ind) # shuffle it
    }
    best_k <- search_window_ind[1]
    best_distance <- 1E10
    
    
    if (dist_metric == "d1"){
      for (k in search_window_ind) {
        candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                  candptind = k,
                                                  lagvectors = lag_vectors,
                                                  ccdiff = cdiff,
                                                  lenout = len_out)
        if(isTRUE(weighted)){
          dist_data_events <- data_event_dist1(cand_fdata = tdata[candidate_ngbs_ind, ],
                                               pfdata = simdata[nb_pts, ], 
                                               wts = 1/lag_vector_norms)
        }else{
          dist_data_events <- data_event_dist1(cand_fdata = tdata[candidate_ngbs_ind, ],
                                               pfdata = simdata[nb_pts, ], 
                                               wts = NULL)
        }
        
        if (dist_data_events < best_distance) {
          best_k <- k
          best_distance <- dist_data_events
        }
      }
    } else if(dist_metric == "d3"){
      for (k in search_window_ind) {
        candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                  candptind = k,
                                                  lagvectors = lag_vectors,
                                                  ccdiff = cdiff,
                                                  lenout = len_out)
        if(isTRUE(weighted)){
          dist_data_events <- data_event_dist_deriv(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                    pfdata = simdata[nb_pts, ],
                                                    wts = 1/lag_vector_norms,
                                                    deriv_order = 1)
        }else{
          dist_data_events <- data_event_dist_deriv(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                    pfdata = simdata[nb_pts, ],
                                                    wts = NULL,
                                                    deriv_order = 1)
        }
        if (dist_data_events < best_distance) {
          best_k <- k
          best_distance <- dist_data_events
        }
      }
    } else if (dist_metric == "d2"){
      for (k in search_window_ind) {
        candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                  candptind = k,
                                                  lagvectors = lag_vectors,
                                                  ccdiff = cdiff,
                                                  lenout = len_out)
        
        dist_data_events <- data_event_dist2(cand_fdata = tdata[candidate_ngbs_ind, ],
                                             pfdata = simdata[nb_pts, ],
                                             dmax = dmax)
        
        if (dist_data_events < best_distance) {
          best_k <- k
          best_distance <- dist_data_events
        }
      }
    } else if (dist_metric == "d4"){
      muod_indices <- matrix(nrow = lswi, ncol = 3)
      if(isTRUE(weighted)){
        for (k in seq_along(search_window_ind)) {
          candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                    candptind = search_window_ind[k],
                                                    lagvectors = lag_vectors,
                                                    ccdiff = cdiff,
                                                    lenout = len_out)
          
          
          
          muod_indices[k, ] <- data_event_sim_map_muod(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                       pfdata = simdata[nb_pts, ], 
                                                       wts = 1/lag_vector_norms)
          
        }
      }else{
        for (k in seq_along(search_window_ind)) {
          candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                    candptind = search_window_ind[k],
                                                    lagvectors = lag_vectors,
                                                    ccdiff = cdiff,
                                                    lenout = len_out)
          muod_indices[k, ] <- data_event_sim_map_muod(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                       pfdata = simdata[nb_pts, ])
          
        }
      }
      muod_indices <- muod_indices/rep(matrixStats::colMaxs(muod_indices), each = lswi)
      if(is.null(wts_muod)){
        muod_indices <- rowMeans(muod_indices)
      }else{
        wts_muod <- wts_muod/sum(wts_muod)
        muod_indices <- rowSums(muod_indices*rep(wts_muod, each = lswi))
      }
      best_k <- search_window_ind[which.min(muod_indices)]
      best_distance <- min(muod_indices)
    } else {
      stop("Error: method unknown \n")
    }
    
    #cat("::--best candidate for location ", i, " is: ", best_k, "::\n")
    # copy the value of the best candidate to simulation grid
    simdata[i, ] <- tdata[best_k,]
    counter <- counter + 1
    cond_points <- c(cond_points,i)
    #cat("::---found a candidate for location: ", i, "-:: \n")
  }
  
  return(simdata)
  
}

##### functional mps with gap filling (with support for inner-outward filling order) #####
# simgrid: coords of the training image and simulation grid (assumed to be the same)
# tdata: training image with data
# simdata: simulation_grid with some points missing to fill
# missing_points : indices of points to fill in simdata
# nnghbs: number of neighbours
# xseq: range vector of x and y locations (square grids considered)
# dist_metric: the distance metric
# weighted: whether to weight disance metric
# wts_muod: weights for distance 4
fun_mps3 <- function(simgrid,
                     tdata,
                     simdata,
                     missing_points,
                     nnghbs, 
                     xseq, 
                     dist_metric, 
                     weighted = T,
                     wts_muod = c(1/3, 1/3, 1/3) # amp, mag, shape,
){
  cdiff <- xseq[2]-xseq[1]
  len_out <- length(xseq)
  
  dm <- dim(tdata);
  nr <- dm[1];
  nc <- dm[2]
  
  # define random path across simulation grid
  rowss <- 1:nr
  cond_points <- rowss[!(rowss %in%  missing_points)]
  
  simpoints <- missing_points
  
  lsmp <- length(simpoints)
  
  #find dmax if distance is d2
  if(dist_metric == "d2"){
    tdata_norm <- sqrt(rowSums(tdata^2))
    tnorm_max <- which.max(tdata_norm)
    tnorm_min <- which.min(tdata_norm)
    dmax <- sqrt(sum((tdata[tnorm_max,] - tdata[tnorm_min,])^2))
  }
  
  # start simulation
  counter <- 1
  
  for (i in simpoints) {
    cat("::---Candidate ", counter, "out of ", lsmp, "for location", i, "-::", "\n")
    coords_i <- simgrid[i, ]
    lcond <- length(cond_points)
    
    # find n nearest neighbours to point i
    coords_cond_points <- simgrid[cond_points, ]
    dist_i_cond_points <- apply(coords_cond_points, 1, euc.dist, x2 = coords_i)
    sorted_dists <- sort(dist_i_cond_points, index.return = T)
    ind_nb_pts <- sorted_dists$ix[1:nnghbs] #
    nb_pts <-   cond_points[ind_nb_pts]
    nb_pts_coords <- simgrid[nb_pts, ]
    
    # define search window in TI
    lag_vectors <- nb_pts_coords - rep(coords_i, each = nnghbs)
    if(isTRUE(weighted)) lag_vector_norms <- sqrt(rowSums(lag_vectors^2))
    xy_max_limit <- round(max(xseq) - matrixStats::colMaxs(lag_vectors), 2)
    xy_min_limit <- round(min(xseq) - matrixStats::colMins(lag_vectors),2)
    test_cond <- (simgrid[, 1] <= xy_max_limit[1]) & # problem of rounding
      (simgrid[, 1] >= xy_min_limit[1]) &
      (simgrid[, 2] <= xy_max_limit[2]) &
      (simgrid[, 2] >= xy_min_limit[2])
    search_window_ind <- which(test_cond) # indexes simgrid
    lswi <- length(search_window_ind)
    if(lswi > 1){ # sampling behaves differently when it receives an integer
      search_window_ind <- sample(search_window_ind) # shuffle it
    }
    best_k <- search_window_ind[1]
    best_distance <- 1E10
    
    
    if (dist_metric == "d1"){
      for (k in search_window_ind) {
        candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                  candptind = k,
                                                  lagvectors = lag_vectors,
                                                  ccdiff = cdiff,
                                                  lenout = len_out)
        if(isTRUE(weighted)){
          dist_data_events <- data_event_dist1(cand_fdata = tdata[candidate_ngbs_ind, ],
                                               pfdata = simdata[nb_pts, ], 
                                               wts = 1/lag_vector_norms)
        }else{
          dist_data_events <- data_event_dist1(cand_fdata = tdata[candidate_ngbs_ind, ],
                                               pfdata = simdata[nb_pts, ], 
                                               wts = NULL)
        }
        
        if (dist_data_events < best_distance) {
          best_k <- k
          best_distance <- dist_data_events
        }
      }
    } else if(dist_metric == "d3"){
      for (k in search_window_ind) {
        candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                  candptind = k,
                                                  lagvectors = lag_vectors,
                                                  ccdiff = cdiff,
                                                  lenout = len_out)
        if(isTRUE(weighted)){
          dist_data_events <- data_event_dist_deriv(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                    pfdata = simdata[nb_pts, ],
                                                    wts = 1/lag_vector_norms,
                                                    deriv_order = 1)
        }else{
          dist_data_events <- data_event_dist_deriv(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                    pfdata = simdata[nb_pts, ],
                                                    wts = NULL,
                                                    deriv_order = 1)
        }
        if (dist_data_events < best_distance) {
          best_k <- k
          best_distance <- dist_data_events
        }
      }
    } else if (dist_metric == "d2"){
      for (k in search_window_ind) {
        candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                  candptind = k,
                                                  lagvectors = lag_vectors,
                                                  ccdiff = cdiff,
                                                  lenout = len_out)
        
        dist_data_events <- data_event_dist2(cand_fdata = tdata[candidate_ngbs_ind, ],
                                             pfdata = simdata[nb_pts, ],
                                             dmax = dmax)
        
        if (dist_data_events < best_distance) {
          best_k <- k
          best_distance <- dist_data_events
        }
      }
    } else if (dist_metric == "d4"){
      muod_indices <- matrix(nrow = lswi, ncol = 3)
      if(isTRUE(weighted)){
        for (k in seq_along(search_window_ind)) {
          candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                    candptind = search_window_ind[k],
                                                    lagvectors = lag_vectors,
                                                    ccdiff = cdiff,
                                                    lenout = len_out)
          
          
          
          muod_indices[k, ] <- data_event_sim_map_muod(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                       pfdata = simdata[nb_pts, ], 
                                                       wts = 1/lag_vector_norms)
          
        }
      }else{
        for (k in seq_along(search_window_ind)) {
          candidate_ngbs_ind <- find_cand_ngbs_inds(simgridd = simgrid,
                                                    candptind = search_window_ind[k],
                                                    lagvectors = lag_vectors,
                                                    ccdiff = cdiff,
                                                    lenout = len_out)
          muod_indices[k, ] <- data_event_sim_map_muod(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                       pfdata = simdata[nb_pts, ])
          
        }
      }
      muod_indices <- muod_indices/rep(matrixStats::colMaxs(muod_indices), each = lswi)
      if(is.null(wts_muod)){
        muod_indices <- rowMeans(muod_indices)
      }else{
        wts_muod <- wts_muod/sum(wts_muod)
        muod_indices <- rowSums(muod_indices*rep(wts_muod, each = lswi))
      }
      best_k <- search_window_ind[which.min(muod_indices)]
      best_distance <- min(muod_indices)
    } else {
      stop("Error: method unknown \n")
    }
    
    #cat("::--best candidate for location ", i, " is: ", best_k, "::\n")
    # copy the value of the best candidate to simulation grid
    simdata[i, ] <- tdata[best_k,]
    counter <- counter + 1
    cond_points <- c(cond_points,i)
    #cat("::---found a candidate for location: ", i, "-:: \n")
  }
  
  return(simdata)
  
}



#### simulation ####

## run fmps simulation (to fill gaps in an SG using a TI),
## for different distances, no of neighbours, and no of conditioning data,
## and number of repetitions
run_simulation2 <- function(no_of_basis = 10, 
                           field_params = list(sigma = rep(1, 10),
                                               beta = rep(.063, 10), 
                                               nu =   rep(1.5, 10)),
                           distance_metric = c("d1", "d2", "d3",
                                               "d4"),
                           no_of_ngbs = c(20, 50, 100),
                           xlocs = round(seq(from = 0, to = 1, len = 51),2), 
                           ylocs = round(seq(from = 0, to = 1, len = 51),2), 
                           weighted = NULL,
                           wts_muod = NULL,
                           npoints_missing = 506,
                           missing_mode = c("box", "scattered"),
                           rep = 5){
  missing_mode <- match.arg(missing_mode)

  # simulate training image if it doesn't already exist
  if(!file.exists("sim_intermediate_data/sgti.RData")){
    ti <- simulate_functional_rf(xseq = xlocs,
                     yseq = ylocs,
                     nobasis = no_of_basis,
                     field_params = field_params)
    # simulate test sim_grid
    sg <- simulate_functional_rf(xseq = xlocs,
                     yseq = ylocs,
                     nobasis = no_of_basis,
                     field_params = field_params)
    save(sg,ti, file = paste0("sim_intermediate_data/sgti.RData"))
  } else{
    load("sim_intermediate_data/sgti.RData")
  }
  # get indices of points to delete
  if(missing_mode == "scattered"){
    set.seed(5)
    missing_indices <- sample.int(nrow(sg$sim_data), size = npoints_missing)
    simdatares <- sg$sim_data
    simdatares[missing_indices, ] <- NA
  }else{
    testcond <- (sg$coord[,1] >= xlocs[1]) &
      (sg$coord[,1] <= xlocs[23]) & 
      (sg$coord[,2] >= ylocs[1]) & 
      (sg$coord[,2] <= ylocs[23])
    missing_indices <- which(testcond)
    simdatares <- sg$sim_data
    simdatares[missing_indices, ] <- NA
  }

  results_list <- list()
  
  for (di in distance_metric) {
      for (nng in no_of_ngbs) {
        for (i in 1:rep) {
          cat("::-Nngbs: ", nng, "; dist: ",di, " rep: ",i, "-::\n")
          simimage <- fun_mps2(simgrid = ti$coords,
                              tdata = ti$sim_data,
                              simdata = simdatares,
                              missing_points = missing_indices,
                              nnghbs = nng,
                              xseq = xlocs,
                              dist_metric = di, 
                              weighted = weighted, 
                              wts_muod = wts_muod)
          intmlist <- list(no_neighbour = nng,
                           distmetric = di,
                           simmimage = simimage,
                           rep = i)

          results_list[[length(results_list)+1]] <- intmlist
          # save(results_list,
          #      file = paste0("sim_intm_data/results_list_",
          #                    "_dist_", di,
          #                    "_nng_", nng,
          #                    "_missing_", missing_mode,
          #                    ".RData"))
        }
        
      }
  }
  final_res <- list(training = ti,
                    simtest = sg,
                    missing_points = missing_indices,
                    results = results_list)
  save(final_res, file = paste0("sim_final_results/final_res_",	
                                "_dist_", distance_metric,
                                "_nng_", no_of_ngbs,
                                "_missing_", missing_mode,
                                ".RData"))
  return(final_res)
  
}

# for inner filling order

fill_along_y <- function(prev_data){
  new_data <- prev_data[-1,]
  new_data$V2 <- new_data$V1
  new_data$V1 <- prev_data$V2[-1]
  return(new_data)
}

fill_along_x <- function(prev_data){
  new_data <- prev_data
  new_data$V1 <- new_data$V2
  new_data$V2 <- new_data$V2[1]
  return(new_data)
}
## run fmps simulation (to fill gaps in an SG using a TI, but with inner_outward_filling),
## for different distances, no of neighbours, and no of conditioning data,
## and number of repetitions
run_simulation3 <- function(no_of_basis = 10, 
                            field_params = list(sigma = rep(1, 10),
                                                beta = rep(.063, 10), 
                                                nu =   rep(1.5, 10)),
                            distance_metric = c("d1", "d2", "d3",
                                                "d4"),
                            no_of_ngbs = c(20, 50, 100),
                            xlocs = round(seq(from = 0, to = 1, len = 51),2), 
                            ylocs = round(seq(from = 0, to = 1, len = 51),2), 
                            weighted = NULL,
                            wts_muod = NULL,
                            npoints_missing = 506,
                            inner_outward = T,
                            missing_mode = c("box", "scattered"),
                            rep = 5){
  missing_mode <- match.arg(missing_mode)
  
  # simulate training image
  if(!file.exists("sim_intermediate_data/sgti.RData")){
    ti <- sim_image2(xseq = xlocs,
                     yseq = ylocs,
                     nobasis = no_of_basis,
                     field_params = field_params)
    # simulate test sim_grid
    sg <- sim_image2(xseq = xlocs,
                     yseq = ylocs,
                     nobasis = no_of_basis,
                     field_params = field_params)
    save(sg,ti, file = paste0("sim_intermediate_data/sgti.RData"))
  } else{
    load("sim_intermediate_data/sgti.RData")
  }
  # get indices of points to delete
  if(missing_mode == "scattered"){
    set.seed(5)
    missing_indices <- sample.int(nrow(sg$sim_data), size = npoints_missing)
    simdatares <- sg$sim_data
    simdatares[missing_indices, ] <- NA
  }else{
    testcond <- (sg$coord[,1] >= xlocs[1]) &
      (sg$coord[,1] <= xlocs[23]) & 
      (sg$coord[,2] >= ylocs[1]) & 
      (sg$coord[,2] <= ylocs[23])
    missing_indices <- which(testcond)
    simdatares <- sg$sim_data
    simdatares[missing_indices, ] <- NA
    if(inner_outward){
      # prep initial data
      ttdt <- as.data.frame(sg$coords[missing_indices, ])
      ds1 <- ttdt[with(ttdt, order(V2, V1, decreasing = T)),]
      ds1 <- ds1[1:which(ds1$V1 == 0)[1],]
      track_data <-comb_data <- ds1
      # rearrange coords inner to outward
      while (nrow(track_data)>1) {
        latest_data <- fill_along_y(track_data)
        track_data <- fill_along_x(latest_data)
        comb_data <- rbind(comb_data, latest_data, track_data)
        cat("::current row number is", nrow(track_data), "\n")
      }
      # convert comb_data to indices in original indices
      arranged_indices <- rep(0,nrow(ttdt))
      for(i in 1:nrow(comb_data)){
        crow <- matrix(rep(unlist(comb_data[i, ]), nrow(ttdt)), byrow = T, nc = 2)
        arranged_indices[i] <-  which(rowSums(ttdt == crow) == 2)
      }
      missing_indices <- missing_indices[arranged_indices]
    }
  }
  # reorder missing points for inner-outward filling
  
  # save(missing_indices, file = paste0("sim_intm_data/gap",
  #                                     "_missing_", missing_mode,
  #                                     "_dist_", distance_metric,
  #                                     "_nng_", no_of_ngbs,
  #                                     ".RData"))
  results_list <- list()
  
  for (di in distance_metric) {
    for (nng in no_of_ngbs) {
      for (i in 1:rep) {
        cat("::-Nngbs: ", nng, "; dist: ",di, " rep: ",i, "-::\n")
        simimage <- fun_mps3(simgrid = ti$coords,
                             tdata = ti$sim_data,
                             simdata = simdatares,
                             missing_points = missing_indices,
                             nnghbs = nng,
                             xseq = xlocs,
                             dist_metric = di, 
                             weighted = weighted, 
                             wts_muod = wts_muod)
        intmlist <- list(no_neighbour = nng,
                         distmetric = di,
                         simmimage = simimage,
                         rep = i)
        
        results_list[[length(results_list)+1]] <- intmlist
        # save(results_list,
        #      file = paste0("sim_intm_data/results_list_",
        #                    "_dist_", di,
        #                    "_nng_", nng,
        #                    "_missing_", missing_mode,
        #                    ".RData"))
      }
      
    }
  }
  final_res <- list(training = ti,
                    simtest = sg,
                    missing_points = missing_indices,
                    results = results_list)
  save(final_res, file = paste0("sim_final_results/final_res_",	
                                "_dist_", distance_metric,
                                "_nng_", no_of_ngbs,
                                "_missing_", missing_mode,
                                "inner_outward",
                                ".RData"))
  return(final_res)
}


#### Distance 1 ####

## box

tt_5 <- run_simulation2(distance_metric = "d1", 
                       no_of_ngbs = 5, 
                       rep = 1, 
                       missing_mode = "box")

# tt_10 <- run_simulation2(distance_metric = "d1",
#                         no_of_ngbs = 10, 
#                         rep = 4, 
#                         missing_mode = "box")
# 
# tt_20 <- run_simulation2(distance_metric = "d1", 
#                         no_of_ngbs = 20, 
#                         rep = 4, 
#                         missing_mode = "box")
# 
# tt_50 <- run_simulation2(distance_metric = "d1", 
#                          no_of_ngbs = 50, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_100 <- run_simulation2(distance_metric = "d1", 
#                           no_of_ngbs = 100, 
#                           rep = 4, 
#                           missing_mode = "box")
#### Distance 2 #####

## box

# tt_5 <- run_simulation2(distance_metric = "d2", 
#                         no_of_ngbs = 5, 
#                         rep = 4, 
#                         missing_mode = "box")
# 
# tt_10 <- run_simulation2(distance_metric = "d2",
#                          no_of_ngbs = 10, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_20 <- run_simulation2(distance_metric = "d2", 
#                          no_of_ngbs = 20, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_50 <- run_simulation2(distance_metric = "d2", 
#                          no_of_ngbs = 50, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_100 <- run_simulation2(distance_metric = "d2", 
#                           no_of_ngbs = 100, 
#                           rep = 4, 
#                           missing_mode = "box")

#### Distance deriv ######

## box

tt_5 <- run_simulation2(distance_metric = "d3", 
                        no_of_ngbs = 5, 
                        rep = 4, 
                        missing_mode = "box")

# tt_10 <- run_simulation2(distance_metric = "d3",
#                          no_of_ngbs = 10, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_20 <- run_simulation2(distance_metric = "d3", 
#                          no_of_ngbs = 20, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_50 <- run_simulation2(distance_metric = "d3", 
#                          no_of_ngbs = 50, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_100 <- run_simulation2(distance_metric = "d3", 
#                           no_of_ngbs = 100, 
#                           rep = 4, 
#                           missing_mode = "box")

#### Distance 4 ####### 

## box
# tt_5 <- run_simulation2(distance_metric = "d4", 
#                         no_of_ngbs = 5, 
#                         rep = 4, 
#                         missing_mode = "box")
# 
# tt_10 <- run_simulation2(distance_metric = "d4",
#                          no_of_ngbs = 10, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_20 <- run_simulation2(distance_metric = "d4", 
#                          no_of_ngbs = 20, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_50 <- run_simulation2(distance_metric = "d4", 
#                          no_of_ngbs = 50, 
#                          rep = 4, 
#                          missing_mode = "box")
# 
# tt_100 <- run_simulation2(distance_metric = "d4", 
#                           no_of_ngbs = 100, 
#                           rep = 4, 
#                           missing_mode = "box")

#### Inner outward ####
tt_5 <- run_simulation3(distance_metric = "d1",
                        no_of_ngbs = 5,
                        rep = 4,
                        missing_mode = "box")

tt_5 <- run_simulation3(distance_metric = "d3",
                        no_of_ngbs = 5,
                        rep = 1,
                        missing_mode = "box")


##### Example plotting

# for plotting all reps of a distance + SG and TI
plot_rf_gap2 <- function(norm_train, norm_testim_missing, norm_simg){
  plotdttt1 <- data.frame(x = rep(ti$coords[,1], 6),
                          y = rep(ti$coords[,2], 6),
                          z =  c(norm_train,
                                 norm_testim_missing,
                                 unlist(norm_simg)), 
                          field = factor(rep(c("training_image",
                                               "test_image_gap",
                                               paste0("simulation", 1:4) ),
                                             each = 2601 ), 
                                         levels = c("training_image",
                                                    "test_image_gap",
                                                    paste0("simulation", 1:4) ))) 
  ggplot(plotdttt1, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    facet_wrap(~field, nrow = 2, ncol = 3)+
    scale_fill_viridis(discrete = FALSE,
                       option = "turbo",
                       #limits = c(-4.5, 4.5),
                       #breaks = c(-4,-3, -2, -1, 0, 1, 2, 3,4),
                       guide = guide_colorbar(label.position = "left",
                                              label.hjust = 1, 
                                              barwidth = 0.5,
                                              barheight = 10, 
                                              draw.ulim = T,
                                              draw.llim = T)) +
    theme_bw()
}

# plot for distance 3

load("sim_final_results/final_res__dist_d3_nng_5_missing_box.RData")
load("sim_intermediate_data/sgti.RData")

norm_training_image <- sqrt(rowSums(final_res$training$sim_data^2))
norm_test_image <- sqrt(rowSums(final_res$simtest$sim_data^2))
norm_test_image_missing <- norm_test_image
norm_test_image_missing[final_res$missing_points]  <- NA
norm_sg_d3 <- list()
for (i in 1:4) {
  norm_sg_d3[[i]] <- sqrt(rowSums(final_res$results[[i]]$simmimage^2))
}


pdf("sim_d3_gap_box.pdf", width = 9.5, height = 7)
plot_rf_gap2(norm_train = norm_training_image,
             norm_testim_missing = norm_test_image_missing, 
             norm_simg = norm_sg_d3)
dev.off()




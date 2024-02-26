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
  # assumed that length of cand_ngbhd_inds and cand_ngbhd_inds are the sam
  
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
  # assumed that length of cand_ngbhd_inds and cand_ngbhd_inds are the same
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
  # assumed that length of cand_ngbhd_inds and cand_ngbhd_inds are the same

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

data_event_sim_map_muod <- function(cand_fdata,
                                    pfdata,
                                    wts = NULL){
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
  # convertion of lag vectors and candidate point ind to negbhd inds formula:
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

##### functional mps #####
# simgrid: coords of the training image and simulation grid (assumed to be the same)
# tdata: training image with data
# simdata: another grid from which initical conditioning data is sampled from
# ncond: number of conditioning data
# nnghbs: number of neighbours
# xseq: range vector of x and y locations (square grids considered)
# dist_metric: the distance metric
# weighted: whether to weight disance metric
# wts_muod: weights for distance 4
fun_mps <- function(simgrid,
                    tdata,
                    simdata,
                    ncond,
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
  
  simdatares <- matrix(NA, nrow = nr, ncol = nc)
  
  if(missing(ncond)){ # unconditional simulation
    
  } else{ # conditional simulation
    cond_points <- sample.int(nr, size = ncond)
    simdatares[cond_points, ] <- simdata[cond_points, ]
  }
  
  # define random path across simulation grid
  rowss <- 1:nr
  simpoints <- rowss[!(rowss %in%  cond_points)]
  if(length(simpoints) > 1){
    simpoints <- sample(simpoints) # shuffle simpoint
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
    
    # check if conditioning points is less nngbs
    if(lcond < nnghbs){
      nnghbs_new <- lcond 
    }else{
      nnghbs_new <- nnghbs
    }
    # find n nearest neighbours to point i
    coords_cond_points <- simgrid[cond_points, ]
    dist_i_cond_points <- apply(coords_cond_points, 1, euc.dist, x2 = coords_i)
    sorted_dists <- sort(dist_i_cond_points, index.return = T)
    ind_nb_pts <- sorted_dists$ix[1:nnghbs_new] #
    nb_pts <-   cond_points[ind_nb_pts]
    nb_pts_coords <- simgrid[nb_pts, ]
    
    # define search window in TI
    lag_vectors <- nb_pts_coords - rep(coords_i, each = nnghbs_new)
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
                                               pfdata = simdatares[nb_pts, ], 
                                               wts = 1/lag_vector_norms)
        }else{
          dist_data_events <- data_event_dist1(cand_fdata = tdata[candidate_ngbs_ind, ],
                                               pfdata = simdatares[nb_pts, ], 
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
                                                    pfdata = simdatares[nb_pts, ],
                                                    wts = 1/lag_vector_norms,
                                                    deriv_order = 1)
        }else{
          dist_data_events <- data_event_dist_deriv(cand_fdata = tdata[candidate_ngbs_ind, ],
                                                    pfdata = simdatares[nb_pts, ],
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
                                             pfdata = simdatares[nb_pts, ],
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
                                                       pfdata = simdatares[nb_pts, ], 
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
                                                       pfdata = simdatares[nb_pts, ])
          
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
    simdatares[i, ] <- tdata[best_k,]
    counter <- counter + 1
    cond_points <- c(cond_points,i)
    #cat("::---found a candidate for location: ", i, "-:: \n")
  }
  
  return(simdatares)
  
}



#### simulation ####

# simulate training image
field_params = list(sigma = rep(1, 10),
                    beta = rep(.063, 10),
                    nu =   rep(1.5, 10))
xlocs = round(seq(from = 0, to = 1, len = 51),2)
ylocs = xlocs
ti <- simulate_functional_rf(xseq = xlocs,
                 yseq = ylocs,
                 nobasis = 10,
                 field_params = field_params)
# simulate another functional rf from which to choose
# initial conditioning data
sg <- simulate_functional_rf(xseq = xlocs,
                 yseq = ylocs,
                 nobasis = 10,
                 field_params = field_params)
# save both in intermediate data folder for later use
save(sg,ti, file = paste0("sim_intermediate_data/sgti_bw.RData"))

## run fmps simulation (generating copies of an existing FRF using the FRF as TI),
## for different distances, no of neighbours, and no of conditioning data,
## and number of repetitions
run_simulation <- function(no_of_basis = 10, 
                           field_params = list(sigma = rep(1, 10),
                                               beta = rep(.063, 10), 
                                               nu =   rep(1.5, 10)),
                           distance_metric = c("d1", "d2", "d3", "d4"),
                           conditioning_data = c(10L, 20L, 50L), 
                           no_of_ngbs = c(5L, 10L, 20L),
                           xlocs = round(seq(from = 0, to = 1, len = 51),2), 
                           ylocs = round(seq(from = 0, to = 1, len = 51),2), 
                           weighted = NULL,
                           wts_muod = NULL,
                           rep = 5){
  
  # load training image and inicial conditioning data image
  load("sim_intermediate_data/sgti_bw.RData")
  # ti <- sim_image2(xseq = xlocs,
  #                  yseq = ylocs,
  #                  nobasis = no_of_basis,
  #                  field_params = field_params)
  # #plot_fields_image(ti)
  # 
  # # simulate test sim_grid
  # sg <- sim_image2(xseq = xlocs,
  #                  yseq = ylocs,
  #                  nobasis = no_of_basis,
  #                  field_params = field_params)
  # #plot_fields_image(sg) 
  # save(sg,ti, file = paste0("sim_intm_data/sgti.RData"))
  results_list <- list()
  
  for (di in distance_metric) {
    for (nco in conditioning_data) {
      for (nng in no_of_ngbs) {
        if(nng >= nco){
          next
        }
        for (i in 1:rep) {
          cat("::-Ncond: ", nco, "; nngbs: ", nng, "; dist: ",di, " rep: ",i, "-::\n")
          simimage <- fun_mps(simgrid = ti$coords,
                              tdata = ti$sim_data,
                              simdata = sg$sim_data,
                              ncond = nco,
                              nnghbs = nng,
                              xseq = xlocs,
                              dist_metric = di, 
                              weighted = weighted, 
                              wts_muod = wts_muod)
          intmlist <- list(no_neighbour = nng,
                           no_cond = nco,
                           distmetric = di,
                           simmimage = simimage,
                           rep = i)
          
          results_list[[length(results_list)+1]] <- intmlist
          save(results_list,
               file = paste0("sim_intermediate_data/results_list_",
                             "dist_", di,
                             "_ncond_", nco,
                             "_nng_", nng,
                             ".RData"))
        }
        
      }
      
    }
    
  }
  
  final_res <- list(training = ti,
                    simtest = sg,
                    results = results_list)
  save(final_res, file = paste0("sim_final_results//final_res_",
                                "dist_", distance_metric,
                                "_ncond_", conditioning_data,
                                "_nng_", no_of_ngbs,
                                ".RData")) 
  return(final_res)
  
}


#### Distance 1 ####

tt_10_5 <- run_simulation(distance_metric = "d1",
                          conditioning_data = 10, 
                          no_of_ngbs = 5, 
                          rep = 4)

# tt_20_5 <- run_simulation(distance_metric = "d1",
#                           conditioning_data = 20, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_50_5 <- run_simulation(distance_metric = "d1",
#                           conditioning_data = 50, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_20_10 <- run_simulation(distance_metric = "d1",
#                            conditioning_data = 20, 
#                            no_of_ngbs = 10, 
#                            rep = 4)
# 
# tt_50_10 <- run_simulation(distance_metric = "d1",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 10, 
#                            rep = 4)
# 
# tt_50_20 <- run_simulation(distance_metric = "d1",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 20, 
#                            rep = 4)
# 
# #### Distance 2 #####
# tt_10_5 <- run_simulation(distance_metric = "d2",
#                           conditioning_data = 10, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_20_5 <- run_simulation(distance_metric = "d2",
#                           conditioning_data = 20, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_50_5 <- run_simulation(distance_metric = "d2",
#                           conditioning_data = 50, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_20_10 <- run_simulation(distance_metric = "d2",
#                            conditioning_data = 20, 
#                            no_of_ngbs = 10, 
#                            rep = 4)
# 
# tt_50_10 <- run_simulation(distance_metric = "d2",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 10, 
#                            rep = 4)
# 
# tt_50_20 <- run_simulation(distance_metric = "d2",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 20, 
#                            rep = 4)
# 
# #### Distance 3 ######
# 
# tt_10_5  <- run_simulation(distance_metric = "d3",
#                            conditioning_data = 10, 
#                            no_of_ngbs = 5, 
#                            rep = 4)
# 
# tt_20_5 <- run_simulation(distance_metric = "d3",
#                           conditioning_data = 20, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_50_5 <- run_simulation(distance_metric = "d3",
#                           conditioning_data = 50, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_20_10  <- run_simulation(distance_metric = "d3",
#                             conditioning_data = 20, 
#                             no_of_ngbs = 10, 
#                             rep = 4)
# 
# tt_50_10 <- run_simulation(distance_metric = "d3",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 10, 
#                            rep = 4)
# 
# tt_50_20 <- run_simulation(distance_metric = "d3",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 20, 
#                            rep = 4)
# 
# #### Distance 4 #####
# tt_10_5  <- run_simulation(distance_metric = "d4",
#                            conditioning_data = 10, 
#                            no_of_ngbs = 5, 
#                            rep = 4)
# 
# tt_20_5 <- run_simulation(distance_metric = "d4",
#                           conditioning_data = 20, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_50_5 <- run_simulation(distance_metric = "d4",
#                           conditioning_data = 50, 
#                           no_of_ngbs = 5, 
#                           rep = 4)
# 
# tt_20_10<- run_simulation(distance_metric = "d4",
#                           conditioning_data = 20, 
#                           no_of_ngbs = 10, 
#                           rep = 4)
# 
# tt_50_10 <- run_simulation(distance_metric = "d4",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 10, 
#                            rep = 4)
# 
# tt_50_20 <- run_simulation(distance_metric = "d4",
#                            conditioning_data = 50, 
#                            no_of_ngbs = 20, 
#                            rep = 4)


#### Example plotting of results ####

# load training grid (ti) and initial conditioning data grid (sg)
load("sim_intermediate_data/sgti_bw.RData")
# load a result
load("sim_final_results/final_res_dist_d1_ncond_10_nng_5.RData")

plot_fields_image <- function(spatial_fn_field, nr = 2, nc = 5){
  nbasis <- length(spatial_fn_field$rndfields)
  dtplot2 <- data.frame(x = rep(spatial_fn_field$coords[,1], times = nbasis), 
                        y = rep(spatial_fn_field$coords[,2], times = nbasis), 
                        z = unlist(spatial_fn_field$rndfields), 
                        field = factor(paste0("field",
                                              rep(1:nbasis, each = nrow(ti$coords))),
                                       levels = paste0("field",1:nbasis)))
  p2 <- ggplot(dtplot2, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    facet_wrap(~field, nrow = nr, ncol = nc)+
    scale_fill_viridis(discrete = FALSE,
                       option = "turbo",
                       limits = c(-4.5, 4.5),
                       breaks = c(-4,-3, -2, -1, 0, 1, 2, 3,4),
                       guide = guide_colorbar(label.position = "left",
                                              label.hjust = 1, 
                                              barwidth = 0.5,
                                              barheight = 10, 
                                              draw.ulim = T,
                                              draw.llim = T)) +
    theme_bw()
  
  p2
  
}

# plot rfs used for sg
plot_fields_image(sg)
# plot rfs used for ti
plot_fields_image(ti)

# plot n reps of simulated FRF
plot_rfs <- function(norm_simg, color_limits, nreps){
  plotdttt1 <- data.frame(x = rep(final_res$training$coords[,1], nreps),
                          y = rep(final_res$training$coords[,2], nreps), 
                          z =  unlist(norm_simg), 
                          field = factor(rep(paste0("Functional_RF", 1:nreps), each = 2601), 
                                         levels =  paste0("Functional_RF", 1:nreps) ))
  
  ggplot(plotdttt1, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    facet_wrap(~field)+
    scale_fill_viridis(discrete = FALSE,
                       option = "turbo",
                       limits = color_limits,
                       breaks = c(5,10, 15, 20),
                       guide = guide_colorbar(label.position = "left",
                                              label.hjust = 1, 
                                              barwidth = 0.5,
                                              barheight = 10, 
                                              draw.ulim = T,
                                              draw.llim = T)) +
    theme_bw()
}

nreps <- 4 # from the number of reps conducted in simulation

norm_training_image <- sqrt(rowSums(final_res$training$sim_data^2))
norm_cond_image <- sqrt(rowSums(final_res$simtest$sim_data^2))

norm_sg_d1 <- list()
for (i in 1:nreps) {
  norm_sg_d1[[i]] <- sqrt(rowSums(final_res$results[[i]]$simmimage^2))
}

lims <- c(min(min(norm_training_image), min(norm_cond_image)),
          max(max(norm_training_image), max(norm_cond_image)))
pdf("sim_new_copy_d1.pdf", width = 9.5, height = 9)
plot_rfs(norm_simg = norm_sg_d1, color_limits = lims, nreps = nreps)
dev.off()

# used to plot four simulated rfs for each distance + train + cond_data
plot_rfs1 <- function(norm_train, norm_cimage, norm_simg){
  plotdttt1 <- data.frame(x = rep(ti$coords[,1], 6),
                          y = rep(ti$coords[,2], 6), 
                          z =  c(norm_train,
                                 norm_cimage, 
                                 unlist(norm_simg)), 
                          field = factor(rep(c("training_image",
                                               "conditioning_data",
                                               paste0("Functional_RF", 1:4) ),
                                             each = 2601 ),
                                         levels = c("conditioning_data",
                                                    "training_image",
                                                    paste0("Functional_RF", 1:4) ))) 
  ggplot(plotdttt1, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    facet_wrap(~field, nrow = 2, ncol = 3)+
    scale_fill_viridis(discrete = FALSE,
                       option = "turbo",
                       limits = c(min(min(norm_train), min(norm_cimage)),
                                  max(max(norm_train), max(norm_cimage))),
                       breaks = c(5,10, 15, 20),
                       guide = guide_colorbar(label.position = "left",
                                              label.hjust = 1, 
                                              barwidth = 0.5,
                                              barheight = 10, 
                                              draw.ulim = T,
                                              draw.llim = T)) +
    theme_bw()
}
# used to plot training + cond_data
plot_rf_tr_ci <- function(norm_train, norm_cimage){
  plotdt <- data.frame(x = rep(final_res$training$coords[,1], 2),
                       y = rep(final_res$training$coords[,2], 2), 
                       z =  c(norm_train, norm_cimage), 
                       field = factor(rep(c("training_image", "conditioning_data"), each = 2601)))
  ggplot(plotdt, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    facet_wrap(~field, nrow = 1, ncol = 2)+
    scale_fill_viridis(discrete = FALSE,
                       option = "turbo",
                       limits = c(min(min(norm_train), min(norm_cimage)),
                                  max(max(norm_train), max(norm_cimage))),
                       breaks = c(5,10, 15, 20),
                       guide = guide_colorbar(label.position = "left",
                                              label.hjust = 1, 
                                              barwidth = 0.5,
                                              barheight = 10, 
                                              draw.ulim = T,
                                              draw.llim = T)) +
    theme_bw()
  
}

# plot of the two rfs: conditiong data and training image
pdf("training_conditioning1.pdf", width = 9, height = 4.5)
plot_rf_tr_ci(norm_train = norm_training_image, norm_cimage = norm_cond_image)
dev.off()

# plots of simulated rfs for d1 + conditioning and training image data
pdf("sim_new_copy_d1_cond_tr.pdf", width = 9.5, height = 6)
plot_rfs1(norm_simg = norm_sg_d1,
          norm_train = norm_training_image,
          norm_cimage = norm_cond_image)
dev.off()

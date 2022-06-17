# Ungrouping and Decumulation of age-specific count data over an
# irregular 2D grid

# Init ------------------------------------------------------------

library(EigenR)
library(tidyverse)
library(plotly)
library(rayshader)

# Data ------------------------------------------------------------

cnst <- within(list(), {
  # grid for grouped data
  grid = rbind(
    expand.grid(
      a = c(0, 20, 40, 60, 80),
      t = c(1990, 1995, 2000, 2005)
    ),
    expand.grid(
      a = c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80),
      t = 2006:2010
    )
  )
  grid_final_age_width = 10
  grid_final_period_width = 1
  max_age = 100
})

# Functions -------------------------------------------------------

# Check for overlap between two rectangles
# defined by (x_i, nx_i, y_i, ny_i) and
# (x_j, nx_j, y_j, ny_j)
# vectorized over rows of x, y, 
DetectOverlap <- function (x, y, nx, ny, p = 0) {
  
  x2 = x+nx+p
  y2 = y+ny+p
  
  N = length(x)
  
  O <- matrix(0, N, N)
  
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      O[i,j] <-
        y[i]  < y2[j] &
        y2[i] > y[j] &
        x[i]  < x2[j] &
        x2[i] > x[j]
    }
  }
  
  O <- which(O==1, arr.ind = T)
  colnames(O) <- c('i', 'j')
  
  return(O)
  
}

PurgeNegativeElements <- function (x, zero_threshold = 1e-8) {
  a <- x
  while (any(a < 0)) {
    posi <- which(a>0)
    nega <- which(a<0)
    na <- a; na[na>=0] <- 0
    b <- a - na
    nna <- na
    for (i in nega) {
      j <- posi[which.min(abs(i-posi))]
      v <- nna[i]
      nna[j] <- v + nna[j]
      nna[i] <- 0
    }
    a <- b + nna
    #cat(sum(a<0), ' ', sum(a), '\n')
  }
  a[a<zero_threshold] <- 0
  return(a)
}

InferCellDimensions <- function (df, time, age, max_age, final_period_width) {
  # width of age groups
  na <- unlist(lapply(
    split(df, df[[time]]), function (x) c(diff(x[[age]]),
                                          max_age-tail(x[[age]], 1))
  ))
  # width of time periods
  nt <- rep(c(diff(unique(G[[time]])), final_period_width),
            c(table(G[[time]])))
  
  return(list(na = na, nt = nt))
}

PlotMatrix <- function (X) {
  XX <- t(X[nrow(X):1,])
  image(XX, xaxt = 'n', yaxt = 'n')
}

# Read ------------------------------------------------------------
 
# fi <- readRDS('dat/input.rds')
# 
# G_fi <-
#   fi %>%
#   filter(
#     Sex == 'b', !Age %in% c('TOT', 'UNK')
#   ) %>%
#   select(
#     a = Age, t = Date, y = Value
#   ) %>%
#   mutate(
#     a = as.integer(a),
#     date = lubridate::dmy(t),
#     t = as.integer(date - min(date))
#   ) %>%
#   arrange(t, a) %>%
#   filter(t %in% c(400:550))

# Prepare and validate input data ---------------------------------

G <- cnst$grid
#G <- G_fi

# ensure that there are no NA's in input grid
stopifnot(!anyNA(G$a))
stopifnot(!anyNA(G$t))

# order in increasing time, age with age varying faster
G <- G[order(G$t, G$a, decreasing = FALSE),]

# infer age*period width from cells
G_cell_dimensions <-
  InferCellDimensions(G, 't', 'a', cnst$max_age, cnst$grid_final_period_width)

G$na <- G_cell_dimensions$na
G$nt <- G_cell_dimensions$nt

# Check input -----------------------------------------------------

# check if there is any overlap between two age*period regions
stopifnot(
  !(nrow(DetectOverlap(x = G$t, y = G$a, nx = G$nt, ny = G$na, p = 0)) > 0)
)

# Construct composition matrix ------------------------------------

# dimensions of grouped data
n_tG = length(unique(G$t)) # number of time period in grouped data
n_G = nrow(G)              # number of time*age cells in grouped data
i_G = 1:n_G                # index over grouped cells

# dimensions of ungrouped grid
max_aU = max(G$a+G$na-1)  # start of highest age in ungrouped data
max_tU = max(G$t+G$nt-1)  # start of highest period in ungrouped data
min_aU = min(G$a)         # start of lowest age group in ungrouped data
min_tU = min(G$t)         # start of lowest period in ungrouped data
a_U = min_aU:max_aU       # ungrouped ages
t_U = min_tU:max_tU       # ungrouped time points
n_aU = length(a_U)        # number of ungrouped age intervals
n_tU = length(t_U)        # number of ungrouped time intervals
n_U = n_aU*n_tU           # number of ungrouped age x time cells
i_U = 1:n_U               # index over ungrouped cells

# G$it <- rep(unique(G$t)-min(G$t)+1, table(G$t))
# G$j  <- rep(1:n_tG, table(G$it))
# G$ia <- unlist(lapply(split(G, ~it), function (X) X$a-min(X$a)+1))
# # time aggregation
# j_pos_t <- (unique(G$t) - min_tU + 1)
# j_width_t <- unlist(lapply(split(G, ~t), function (x) x[1,'nt']))
# # time composition
# C_t <- matrix(0, n_tU, n_tG)
# for (j in 1:n_tG) {
#   #start <- j_pos_t[j]
#   # if cumulative along time
#   start <- 1
#   end <- j_pos_t[j]+j_width_t[j]-1
#   C_t[start:end, j] <- 1
# }
# # age aggregation
# C <- NULL
# for (j in 1:n_tG) {
#   X <- G[G[,'j']==j,]
#   C_a <- matrix(0, nrow(X), n_aU)
#   for (i in 1:nrow(X)) {
#     start <- X[i,][['ia']]
#     end <- start+X[i,][['na']]-1
#     C_a[i, start:end] <- 1
#   }
#   C <- rbind(C, kronecker(t(C_t[,j]), C_a))
# }

# create empty composition matrix
# number of rows equal to number of grouped time*age cells
# number of columns equal to number of ungrouped time*age cells
C <- matrix(0, nrow = n_G, ncol = n_U)
# column indices of first column where C is 1 for each row
j_start <- (G[['t']]-min_tU)*n_aU + G[['a']] - min_aU + 1
# transform to 1D index (row major order!)
indx <- j_start+(i_G-1)*n_U
# those repeat and get shifted according to the number of fine grid
# age groups with every additional year in the bin
k <- unlist(sapply(G[['nt']], FUN = function (x) 1:x))-1
indx <- rep.int(indx, times = G[['nt']]) + k*n_aU
# cumulate across time
indx <- unlist(sapply(indx, function (i) seq(i, floor(i/n_U)*n_U+1, -n_aU)))
# every starting j gets extended to the right by the corresponding
# width of the age group - 1
tmp_i <- ceiling(indx/n_U) # i position
k <- unlist(sapply(G[['na']][tmp_i], FUN = function (x) 1:x))-1
indx <- rep.int(indx, times = G[['na']][tmp_i]) + k
# change to column major order index
tmp_i <- ceiling(indx/n_U) # i position
tmp_j <- indx-(tmp_i-1)*n_U  # j position
indx_colmaj <- (tmp_j-1)*n_G + tmp_i

C[indx_colmaj] <- 1

# plot the composition matrix
PlotMatrix(C)

# Simulate data ---------------------------------------------------

# quadratic Poisson surface
sim <- list()
sim$grid <- expand.grid(
  seq(-1, 1, length.out = n_aU),
  seq(-1, 1, length.out = n_tU)
)
sim$grid$lambda <- exp(-(sim$grid$Var1^2 + sim$grid$Var2^2))*1e2
sim$grid$y <- rpois(n_U, sim$grid$lambda)
sim$Y <-
  matrix(sim$grid$y, nrow = n_aU, ncol = n_tU, dimnames = list(a_U, t_U))
sim$Lambda <- matrix(sim$grid$lambda, nrow = n_aU, ncol = n_tU,
                     byrow = FALSE)

# the observed counts over the grouped grid
y <- unlist(C %*% c(sim$Y))
G[['y']] <- as.vector(y)

# Actual data -----------------------------------------------------

y <- G[['y']]

# plot grouped data grid
G %>%
  mutate(t = t+nt/2, a = a+na/2) %>%
  ggplot() +
  geom_tile(aes(x = t, y = a, width = nt, height = na, fill = y),
            color = 'black', size = 0.1) +
  scale_fill_viridis_c()

# Setup 2D B-spline basis -----------------------------------------

library(splines)

B_a <- bs(a_U, df = ceiling(n_aU/30),
          degree = 2, intercept = TRUE)
B_t <- bs(t_U, df = ceiling(n_tU/30),
          degree = 2, intercept = TRUE)

# combine 1D spline bases into 2D basis
B_ta <- kronecker(B_t, B_a)

PlotMatrix(B_a)
PlotMatrix(B_t)
PlotMatrix(B_ta)

# Setup linear predictor ------------------------------------------

XX <- NULL

# Complete design matrix ------------------------------------------

X <- cbind(B_ta, XX)

# Moore-Penrose approach ------------------------------------------

# Moore-Penrose inverse of composition matrix...
C_inv <- EigenR::Eigen_pinverse(C)
# ...used to infer the ungrouped and decumulated number of counts
# under the assumption that counts are uniformly distributed
# over each grouped region
y_u_mp <- C_inv %*% y
# redistribute negative part of negative elements to neighbors
# so that negative elements become 0 but the vector sum remains
# identical
y_u_mp <- PurgeNegativeElements(y_u_mp)
# because we're dealing with cumulated counts, the data may be
# left censored, i.e. the starting counts of the cumulative series
# are only known to have occured before the starting point
# we treat those starting counts as an offset
y_a0_mp <- y_u_mp[1:n_aU]
y_a0_mp_offset <- c(y_a0_mp, rep(0, (n_tU-1)*n_aU))
y_u_mp[1:n_aU] <- 0
Y_u_mp <- matrix(y_u_mp, n_aU, n_tU)

PlotMatrix(Y_u_mp)
rgl::persp3d(
  x = a_U, y = t_U-min(t_U), Y_u_mp,
  color = 'grey'
)

# Setup penalty ---------------------------------------------------

# smoothing parameter
lambda <- 1e-6
# difference matrix
# penalize second differences
# between neighboring ages within the same year
# via the corresponding spline coefficients
D_a <- kronecker(diag(ncol(B_t)), diff(diag(ncol(B_a)), diff = 2))
# penalize second differences
# between neighboring years within the same age
D_t <- kronecker(diag(ncol(B_a)), diff(diag(ncol(B_t)), diff = 2))
D <- rbind(D_a, D_t)

# Estimate model coefficients -------------------------------------

ll <- 1
b <- coef(lm(log(y_u_mp+1e-8)~-1+X))

# Perform the iterations
for (i in 1:30) {
  
  b_last <- b       # b from last iteration
  ll_last <- ll   # deviance from last iteration
  
  eta <- X %*% b    # linear predictor
  gam <- exp(eta)+y_a0_mp_offset   # ungrouped mean prediction
  mu <- C %*% gam   # grouped mean prediction    
  ll <- -sum(y*log(y/(mu+1e-8)+1e-8)-(y-mu)) + lambda*sum(D%*%b)
  
  if (is.nan(ll)) {
    while (is.nan(ll)) {
      b2 <- (b_last+b)/2
      b_last <- b; b <- b2
      
      eta <- X %*% b    # linear predictor
      gam <- exp(eta)+y_a0_mp_offset   # ungrouped mean prediction
      mu <- C %*% gam   # grouped mean prediction    
      ll <- -sum(y*log(y/mu+1e-8)-(y-mu)) + lambda*sum(D%*%b)
      
    }
  }
  
  Gam <- gam %*% rep(1, ncol(X)) # include offset
  Q <- C %*% (Gam * X)
  # pseudo observations with added 0's for the penalties, e.g.
  # Db=0 seeks to minimize the differences between the beta's
  z <- c(y - mu + Q %*% b, rep(0, nrow(D)))
  
  w <-        # weights for Poisson likelihood and difference penalty
    c(1/mu, rep(lambda, nrow(D)))
  # estimate new beta's under the penalty
  Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F) 
  b <- Fit$coef

  
  db <- abs(ll - ll_last)
  cat(db, '\n')
  if (db < 1e-4) break
  
    
}

# Predict from model ----------------------------------------------

eta <- X%*%b
y_u_pclm <- exp(eta)
Y_u_pclm <- matrix(y_u_pclm, nrow = n_aU, ncol = n_tU)

# Prepare for export ----------------------------------------------

method <- 'ginv'

# Moore-Penrose approach
if (identical(method, 'ginv')) {
  ungrouped_counts <- list(
    y = y_u_mp + y_a0_mp_offset,
    Y = Y_u_mp + y_a0_mp_offset,
    y_no_offset = y_u_mp,
    y_no_offset = Y_u_mp,
    offset = y_a0_mp_offset
  )
}
# Penalized composite link model
if (identical(method, 'pclm')) {
  ungrouped_counts <- list(
    y = y_u_pclm,
    Y = Y_u_pclm + y_a0_mp_offset,
    y_no_offset = y_u_pclm,
    y_no_offset = Y_u_pclm,
    offset = y_a0_mp_offset
  )
}

# Visualize ungrouped counts --------------------------------------

y_u_df <- expand.grid(
  a = a_U,
  t = t_U
)

y_u_df[['y']] <- ungrouped_counts[['y']]

y_u_df %>%
  mutate(t = t+0.5, a = a+0.5) %>%
  ggplot() +
  geom_tile(aes(x = t, y = a, width = 1, height = 1, fill = y),
            color = 'black', size = 0.1) +
  scale_fill_viridis_c(trans = 'log1p')

# Visualize ungrouped cumulative counts ---------------------------

y_u_df[['y_cumt']] <- c(t(apply(
  ungrouped_counts[['Y']]+ungrouped_counts[['offset']], 1, cumsum
)))

y_u_df <-
  y_u_df %>%
  mutate(a = cut(a, seq(0, 100, 5), right = FALSE)) %>%
  group_by(t, a) %>%
  summarise(y = sum(y_cumt))

y_u_df %>%
  ggplot(aes(x = t, y = y)) +
  geom_step() +
  facet_wrap(~a, scales = 'free_y')

# Calculate re-aggregation error ----------------------------------

y_reagg <- C%*%ungrouped_counts[['y']]

plot(y-y_reagg)
mean(abs(y-y_reagg))

# error by time*age cell
G %>%
  mutate(y2 = y_reagg) %>%
  group_by(t) %>%
  summarise(y = sum(y), y2 = sum(y2), err = abs(y2-y)/y*100) %>%
  arrange(-err)
# total count error
abs(sum(y)-sum(y_reagg))/sum(y)*100

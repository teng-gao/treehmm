############ allele HMM ############

#' Beta-binomial distribution density function
#' A distribution is beta-binomial if p, the probability of success, 
#' in a binomial distribution has a beta distribution with shape 
#' parameters α > 0 and β > 0
#' For more details, see extraDistr::dbbinom
#'
#' @param x vector of quantiles
#' @param size number of trials (zero or more)
#' @param alpha numeric (default=1)
#' @param beta numeric (default=1)
#' @param log boolean (default=FALSE)
#' @return density values returned as numeric vector
#' @examples
#' xx <- 1:1000
#' dbbinom(xx, 1000, 5, 13)
#'
#' @export
dbbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
    cppdbbinom(x, size, alpha, beta, log[1L])
}

#' @keywords internal  
makedensity <- function(distn){
    ## https://github.com/cran/HiddenMarkov/blob/master/R/makedensity.R
    dname <- paste("d", distn[1], sep="")
    x <- paste("function(x, pm, pn=NULL, log=FALSE)
         do.call(\"", dname, "\", c(list(x=x), pm, pn,", sep="")
    if (distn[1]=="glm") x <- paste(x, " list(family=\"", distn[2],
         "\", link=\"", distn[3], "\"),", sep="")
    eval(parse(text=paste(x, " list(log=log)))", sep="")))
}


#' @keywords internal
getj <- function(x, j){
    ## https://github.com/cran/HiddenMarkov/blob/master/R/getj.R
    #   get the jth "row" from a list
    if (is.null(x)) return(NULL)
    n <- length(x)
    for (i in 1:n)
        x[[i]] <- x[[i]][j]
    return(x)
}

############ time inhomogenous univariate HMM ############

#' Viterbi algorithm for allele HMM
#' @keywords internal
viterbi_allele <- function (obj, ...){
#     print('Solving univariate nonhomogenous markov chain')
    x <- obj$x
    dfunc <- makedensity(obj$distn)
    n <- length(x)
    m <- nrow(obj$Pi[[1]])
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)

    nu[1, ] <- log(obj$delta) 

    if (!is.na(x[1])) {
        nu[1, ] <- nu[1, ] + dfunc(x=x[1], obj$pm, getj(obj$pn, 1), log=TRUE)
    }

    if (is.null(obj$logPi)) {
        logPi <- lapply(obj$Pi, log)
    } else {
        logPi = obj$logPi
    }

    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        nu[i, ] <- apply(matrixnu + logPi[[i]], 2, max)
        if (!is.na(x[i])) {
            nu[i, ] <- nu[i, ] + dfunc(x=x[i], obj$pm, getj(obj$pn, i), log=TRUE)
        }
    }
#     if (any(nu[n, ] == -Inf)) 
#         stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    # double check this index of logPi
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[[i + 1]][, y[i + 1]] + nu[i, ])

    LL = max(nu[n, ])
        
    return(y)
}


#' Allele treeHMM with one theta state
#' @param pAD
#' @param DP
#' @param p_s
#' @param Q node marginals
#' @param Q_pair edge marginals
#' @param pa parent chain
#' @param ch children chains
#' @keywords internal
get_allele_treehmm = function(pAD, DP, p_s, Q = NULL, Q_pair = NULL, pa = NULL, ch = NULL, t = 1e-5, theta_min = 0.08, gamma = 20, prior = NULL) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('More than one gamma parameter')
    }
    
    # states
    states = c("neu", "theta_up", "theta_down")

    # no parent or children
    if (is.null(pa) & is.null(ch)) {
        calc_trans_mat = function(p_s, t) {
            matrix(
                c(
                    (1-t), t/2, t/2,
                    t, (1-t)*(1-p_s), (1-t)*p_s, 
                    t, (1-t)*p_s, (1-t)*(1-p_s)
                ),
                ncol = 3,
                byrow = TRUE
            )
        }
        As = lapply(p_s, function(p_s) {calc_trans_mat(p_s, t)})
        logAs = As %>% purrr::map(log)
    } else {
        logAs = treehmm_trans_mat(t, p_s, Q, Q_pair, pa, ch)
        As = logAs %>% purrr::map(exp)
    }
    
    # intitial probabilities
    if (is.null(prior)) {
        prior = rep(1/3, 3)
    }

    alpha_up = (0.5 + theta_min) * gamma
    beta_up = (0.5 - theta_min) * gamma
    alpha_down = beta_up
    beta_down = alpha_up
    alpha_neu = gamma/2
    beta_neu = gamma/2
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha=c(alpha_neu, alpha_up, alpha_down), beta=c(beta_neu, beta_up,beta_down)),
        pn = list(size = DP),
        discrete = TRUE)

    hmm$states = states
    hmm$logPi = logAs

    class(hmm) = 'dthmm.inhom'

    return(hmm)
}

# generate lookup table for conditional state probablities
# p_z[t, z_t_pa, z_t-1, z_t]
get_p_z = function(t, p_s) {

    # calculate state conditionals
    get_z_conditional = function(t, p_s, z_pa, z_t_1, z_t) {

        case_when(
            z_pa == 1 & z_t_1 == 1 & z_t == 1 ~ (1-t), 
            z_pa == 1 & z_t_1 == 1 & z_t == 2 ~ t/2, 
            z_pa == 1 & z_t_1 == 1 & z_t == 3 ~ t/2,
            z_pa == 1 & z_t_1 == 2 & z_t == 1 ~ t,
            z_pa == 1 & z_t_1 == 2 & z_t == 2 ~ (1-t)*(1-p_s),
            z_pa == 1 & z_t_1 == 2 & z_t == 3 ~ (1-t)*p_s,
            z_pa == 1 & z_t_1 == 3 & z_t == 1 ~ t, 
            z_pa == 1 & z_t_1 == 3 & z_t == 2 ~ (1-t)*p_s, 
            z_pa == 1 & z_t_1 == 3 & z_t == 3 ~ (1-t)*(1-p_s),
            # z_pa != 1 & z_t_1 == 1 & z_t == 1 ~ (1-t*eta),
            # z_pa != 1 & z_t_1 == 1 & z_t == 2 ~ t*eta/2, 
            # z_pa != 1 & z_t_1 == 1 & z_t == 3 ~ t*eta/2, 
            # z_pa != 1 & z_t_1 == 1 & z_t == 1 ~ t,
            # z_pa != 1 & z_t_1 == 1 & z_t == 2 ~ (1-t)/2, 
            # z_pa != 1 & z_t_1 == 1 & z_t == 3 ~ (1-t)/2, 
            z_pa != 1 & z_t_1 == 1 & z_t == 1 ~ 1/3,
            z_pa != 1 & z_t_1 == 1 & z_t == 2 ~ 1/3, 
            z_pa != 1 & z_t_1 == 1 & z_t == 3 ~ 1/3, 
            z_pa != 1 & z_t_1 == 2 & z_t == 1 ~ t,
            z_pa != 1 & z_t_1 == 2 & z_t == 2 ~ (1-t)*(1-p_s),
            z_pa != 1 & z_t_1 == 2 & z_t == 3 ~ (1-t)*p_s,
            z_pa != 1 & z_t_1 == 3 & z_t == 1 ~ t, 
            z_pa != 1 & z_t_1 == 3 & z_t == 2 ~ (1-t)*p_s, 
            z_pa != 1 & z_t_1 == 3 & z_t == 3 ~ (1-t)*(1-p_s)
        )
    }

    states = c("neu", "theta_up", "theta_down")

    states_grid = expand.grid(1:3, 1:3, 1:3) %>%
        setNames(c('z_t_pa', 'z_t_1', 'z_t'))

    p_z = sapply(
            states_grid %>% split(1:nrow(.)),
            function(Z){
                get_z_conditional(
                    t,
                    p_s = p_s,
                    Z$z_t_pa,
                    Z$z_t_1,
                    Z$z_t)
            }
        ) %>%
        array(dim = c(length(p_s),3,3,3))

    return(p_z)
}

# calculate transition matrices based on new marginals
#' @param pa parent index
#' @param ch children indices
#' @param Q Q[i, 1:t, z_t]
#' @param Q_pair Q_pair[i, t, z_t-1, z_t]
treehmm_trans_mat = function(t, p_s, Q, Q_pair, pa, ch) {

    bigT = length(p_s)

    eta = 1e2

    p_z = get_p_z(t, p_s)

    # defined for t=2:T
    logf_it = function(pa, ch, t, z_t_1, z_t) {

        if (!is.null(pa)) {
            # p_z[t, , z_t_1, z_t] is z_pa x z_t; Q[pa, t, ] is z_pa x 1 => m_pa is 1 x z_t
            m_pa = colSums(log(p_z[t, , z_t_1, z_t]) * Q[pa, t, ])
        } else {
            m_pa = log(p_z[t, 1, z_t_1, z_t])
        }

        if (!is.null(ch)) {
            m_ch = sapply(
                ch,
                function(i) {
                    sapply(
                        1:3,
                        function(w) {
                            sapply(
                                1:3,
                                function(y) {
                                    log(p_z[t, z_t, y, w]) * Q_pair[i, t, y, w]
                                }
                            ) %>% rowSums
                        }
                    ) %>% rowSums
                }
            ) %>%
            rowSums
            
        } else {
            m_ch = 0
        }

        m_pa + m_ch
    }

    # generate a list of transition matrices across time
    logf_z = lapply(
        1:bigT,
        function(t) {
            if (t == 1) {
                return(matrix(rep(NA,9), ncol = 3, nrow = 3))
            } else {
                sapply(1:3, function(z_t_1) {
                    logf_it(
                        pa = pa,
                        ch = ch,
                        t = t,
                        z_t_1 = z_t_1,
                        z_t = 1:3
                    )
                })
            }
        }
    )

    logf_z_old = logf_z

    # renormalize the transition matrix
    log_gz = rep(0,3)

    for (t in bigT:2) {
        logf_z[[t]] = t(t(logf_z[[t]]) + log_gz)
        log_gz = sapply(1:3, function(i){logSumExp(logf_z[[t]][i,])})
        logf_z[[t]] = logf_z[[t]] - log_gz
    }

    return(logf_z)
}

run_treehmm = function(bulks, pa_dict, ch_dict, theta_min = 0.08, gamma = 20, t = 1e-5, max_iter = 10) {
    
    bulks = bulks %>% split(.$sample)
    
    I = length(bulks)
    N = nrow(bulks[[1]])
    
    k = 1
    
    Q = array(rep(NA, I*N*3), dim = c(I, N, 3))
    Q_pair = array(rep(NA, I*N*3*3), dim = c(I, N, 3, 3))
    S = array(rep(NA, max_iter*I*N), dim = c(max_iter, I, N))
    P = array(rep(NA, max_iter*I*N), dim = c(max_iter, I, N))

    logprob = array(rep(NA, I*N*3), dim = c(I, N, 3))
    logtheta = log(get_p_z(t, bulks[[1]]$p_s))
    F = c()
    
    while (k <= max_iter) {

        message(glue('iter {k}'))
        
        ### update variational distribution ###
        for (i in I:1) {

            if (k == 1) {
                pa = NULL
                ch = NULL
            } else {
                pa = pa_dict[[i]]
                ch = ch_dict[[i]]
            }
            
            treehmm = bulks[[i]] %>% {
                get_allele_treehmm(
                    .$pAD, .$DP, .$p_s, 
                    t = t,
                    theta_min = theta_min,
                    Q = Q, 
                    gamma = gamma, 
                    Q_pair = Q_pair, 
                    pa = pa,
                    ch = ch
                    # prior = c(1-2*t, t, t)
                )
            }
            
            fb = forward_back_allele_R(treehmm)

            Q[i,,] = fb$marginals
            Q_pair[i,,,] = fb$edge_marginals
            logprob[i,,] = fb$logprob
            
            S[k,i,] = viterbi_allele(treehmm)
            P[k,i,] = Q[i,,1]
            
        }


        ### calculate free energy ###
        Fk = 0

        for (i in I:1) {

            pa = pa_dict[[i]]

            # entropy term
            H = sapply(
                2:N,
                function(t) {
                    q_cond = Q_pair[i, t, , ]/Q[i, t-1, ]
                    sum(Q_pair[i, t, , ] * log(q_cond))
                }
            ) %>% sum

            H = H + sum(Q[i,1,] * log(Q[i,1,]))

            # likelihood term
            L = sapply(
                2:N,
                function(t) {
                    
                    emission_term = sum(Q[i,t,] * logprob[i,t,])
                    
                    transition_term = sum(sapply(1:3,
                        function(z_pa) {
                            sum(logtheta[t, z_pa, , ] * Q_pair[i, t, , ] * Q[pa, t, z_pa])
                        })
                    )
                    
                    emission_term + transition_term
                }
            ) %>% sum

            # free energy
            Fk = Fk + H - L

        }

        F = c(F, Fk)

        k = k + 1
    }

    return(list(S = S, P = P, F = F))
}



#' @keywords internal
forward_back_allele = function (obj, ...) {

    # case of one-data point
    if (length(obj$x) == 1) {
        return(NA)
    }
    
    x <- obj$x
    p_x <- makedensity(obj$distn)
    
    m <- nrow(obj$Pi[[1]])
    n <- length(x)
    
    logprob = sapply(1:m, function(k) {
        
        l_x = p_x(x = x,  getj(obj$pm, k), obj$pn, log = TRUE)
        
        l_x[is.na(l_x)] = 0
        
        return(l_x)
        
    })
        
    logphi <- log(as.double(obj$delta))

    if (is.null(obj$logPi)) {
        logPi <- lapply(obj$Pi, log)
    } else {
        logPi = obj$logPi
    }
            
    marginals = forward_backward_compute(obj, logphi, logprob, logPi, n, m)

    colnames(marginals) = obj$states

    return(marginals)
}

forward_back_allele_R = function (obj, ...) {

    # case of one-data point
    if (length(obj$x) == 1) {
        return(NA)
    }
    
    x <- obj$x
    p_x <- makedensity(obj$distn)
    
    m <- nrow(obj$Pi[[1]])
    n <- length(x)
    
    logprob = sapply(1:m, function(k) {
        
        l_x = p_x(x = x,  getj(obj$pm, k), obj$pn, log = TRUE)
        
        l_x[is.na(l_x)] = 0
        
        return(l_x)
        
    })
        
    logphi <- log(as.double(obj$delta))
    logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
    lscale <- as.double(0)

    if (is.null(obj$logPi)) {
        logPi <- lapply(obj$Pi, log)
    } else {
        logPi = obj$logPi
    }

    for (t in 1:n) {
        
        if (t > 1) {
            logphi <- sapply(1:m, function(j) matrixStats::logSumExp(logphi + logPi[[t]][,j]))
        }
                          
        logphi <- logphi + logprob[t,]
                          
        logsumphi <- matrixStats::logSumExp(logphi)
                          
        logphi <- logphi - logsumphi
                          
        lscale <- lscale + logsumphi
                          
        logalpha[t,] <- logphi + lscale
                          
        LL <- lscale
    }

    logbeta <- matrix(as.double(rep(0, m * n)), nrow = n)
    logphi <- log(as.double(rep(1/m, m)))
    lscale <- as.double(log(m))

    for (t in seq(n-1, 1, -1)){
        
        logphi = sapply(1:m, function(i) matrixStats::logSumExp(logphi + logprob[t+1,] + logPi[[t+1]][i,]))

        logbeta[t,] <- logphi + lscale

        logsumphi <- matrixStats::logSumExp(logphi)

        logphi <- logphi - logsumphi

        lscale <- lscale + logsumphi
    }

    marginals = exp(logalpha + logbeta - LL)
    colnames(marginals) = obj$states

    edge_marginals = lapply(
        1:(n-1),
        function(t) {
            x = exp(outer(logalpha[t,], logbeta[t+1,], FUN = '+') + logPi[[t+1]] + logprob[t+1,] - LL)
            return(x)
        }
    )
    # defined for (t-1, t) for t=2:T
    edge_marginals = c(matrix(rep(NA,9), ncol = 3, nrow = 3), edge_marginals)

    edge_marginals = aperm(array(unlist(edge_marginals), dim = c(m,m,n), dimnames = list(obj$states, obj$states, 1:(n))), c(3,1,2))

    return(list(marginals = marginals, edge_marginals = edge_marginals, logalpha = logalpha, logbeta = logbeta, logprob = logprob, LL = LL, logPi = logPi))
}


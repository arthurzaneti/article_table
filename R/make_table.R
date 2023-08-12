#___________________________BASIC_FUNCTIONS_____________________________________

rchen_rpr <- function(lambda, mds, tau=0.5){
  nelementos <- length(mds)
  vetor_rchen <- (log(1 - log(1-stats::runif(nelementos))* ((1-exp(mds^lambda))/log(1-tau))))^(1/lambda) #função quantílica da Chen reparametrizada em termos do quantil
  return (vetor_rchen)
}

log_vero_chen <- function(theta, tau, covariables, y){
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)

  log_vero <- suppressWarnings(log(log(1 - tau) / (1 - exp(mds^lambda))) +
                                 log(lambda) + (lambda - 1) * log(y) +
                                 (log(1 - tau) / (1 - exp(mds^lambda))) * (1 - exp(y^lambda)) + (y^lambda))
  return(sum(log_vero))
}

#_______________________________ESTIMATION______________________________________

estim <- function(vetor_random_chen, covariables, tau=0.5, full = F){
  theta_guesses <- c(rep(1, ncol(covariables) +1))
  tryCatch({suppressWarnings(estimacao <- stats::optim(par = theta_guesses,
                                                fn = log_vero_chen,
                                                tau = tau,
                                                covariables=covariables,
                                                y = vetor_random_chen,
                                                method ="BFGS",
                                                hessian = full,
                                                control = list(fnscale=-1)))
    if(full){
      return(estimacao)
    }else{
      return(estimacao$par)
    }
  },error=function(e){
    return(NULL)
  })
}

eval_estim <- function(n_rvalues, monte_carlo_iterations, theta, hist=F){

  covariables <- matrix(rep(1, n_rvalues))
  for(i in 2:(length(theta) - 1)){
    covariables <- cbind(covariables, matrix(stats::runif(n_rvalues)))
  }
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)
  print(mean(mds))
  estims <- matrix(0, nrow= monte_carlo_iterations, ncol=length(theta))
  pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent",
                         total = monte_carlo_iterations)

  errors <- 0
  for (i in 1:monte_carlo_iterations){
    y <- rchen_rpr(lambda, mds)
    par_hats <- suppressWarnings(estim(y, covariables, tau = 0.5))
    if(is.null(par_hats)){
      errors <- errors + 1
      if(i !=1){
        estims[i,] = estims[i-1,]
      }
    }else{
      estims[i,] = par_hats
    }
    pb$tick()
  }
  eval <- matrix(0, nrow=4, ncol=length(theta))
  eval[1,] <- colMeans(estims)
  eval[2,] <- eval[1,] - theta
  eval[3,] <- apply(estims, 2, stats::sd)
  eval[4,] <- eval[2,]^2 + apply(estims, 2, stats::var)

  cols <- c("\u03BB hat")
  for (i in 1:length(betas)){
    cols <- c(cols, paste0("beta", i, " hat"))
  }
  colnames(eval) <- cols
  rownames(eval) <- c("Mean", "Bias", "Standard Error", "MSE")

  if(hist){
    opar <- graphics::par()
    rows <- floor(sqrt(length(theta) * (16/9)))
    cols <- ceiling(length(theta) / rows)
    graphics::par(mfrow = c(rows, cols), family = "serif")
    hist(estims[, 1], main = "\u03BB hat", xlab = "", ylab = "")
    for(i in 2:length(theta)){
      hist(estims[, i], main = paste0("beta", i-1, " hat"), xlab = "", ylab = "")
    }
    suppressWarnings(graphics::par(opar))
  }

  return (eval)
}
#_____________________________CONFIDENCE_INTERVALS______________________________

ci <- function(n_rvalues, theta, alpha, monte_carlo = FALSE, monte_carlo_hist = FALSE) {
  covariables <- matrix(rep(1, n_rvalues))
  for(i in 2:(length(theta) - 1)){
    covariables <- cbind(covariables, matrix(stats::runif(n_rvalues)))
  }
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)

  y <- rchen_rpr(lambda, mds)

  tryCatch({
    estimation <- estim(y, covariables, full = TRUE)

    inf <- solve(-estimation$hessian)
    par <- estimation$par

    lb <- par - stats::qnorm(1 - alpha/2) * sqrt(diag(inf))
    ub <- par + stats::qnorm(1 - alpha/2) * sqrt(diag(inf))

    if (monte_carlo_hist) {
      results <- list(lb[1] <= theta[1] && theta[1] <= ub[1])
      for (i in 2:length(theta)){
        results <- append(results,(lb[i] <= theta[i] && theta[i] <= ub[i]))
      }
      results <- append(results, lb)
      results <- append(results, ub)
      return(results)
    }
    else if (monte_carlo) {
      results <- c(lb[1] <= theta[1] && theta[1] <= ub[1])
      for (i in 2:length(theta)){
        results <- append(results,(lb[i] <= theta[i] && theta[i] <= ub[i]))
      }
      return(results)
    }
    else {
      return(matrix(c(lb,ub), ncol=2))
    }
  }, error = function(e) {
    return(NULL)
  })
}



eval_ci <- function(n_rvalues, monte_carlo_iterations, theta, alpha, hist = F) {
  covariables <- matrix(rep(1, n_rvalues))
  for(i in 2:(length(theta) - 1)){
    covariables <- cbind(covariables, matrix(stats::runif(n_rvalues)))
  }
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)

  TF_table <- matrix(0, nrow= monte_carlo_iterations, ncol=length(theta))
  errors <- 0

  if (!hist) {
    pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", total = monte_carlo_iterations)

    for (i in 1:monte_carlo_iterations) {
      results <- ci(n_rvalues, theta, alpha, monte_carlo = T)

      if (is.null(results)) {
        if (i != 1) {
          errors <- errors + 1
          TF_table[i-1,] <- TF_table[i-1,]
        }
      } else {
        TF_table[i,] <- results
      }
      pb$tick()
    }

  } else {
    lbs <- matrix(0, monte_carlo_iterations, length(theta))
    ubs <- matrix(0, monte_carlo_iterations, length(theta))
    errors <- 0
    pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", total = monte_carlo_iterations)

    for (i in 1:monte_carlo_iterations) {
      results <- ci(n_rvalues, theta, alpha, monte_carlo_hist = T)

      if (is.null(results)) {
        if (i != 1) {
          TF_table[i,] <- TF_table[i-1,]
          lbs[i,] <- lbs[i-1,]
          ubs[i,] <- ubs[i-1,]
        }
      } else {
        s <- length(theta)
        TF_table[i,] <- unlist(results[1:s])
        lbs[i, ] <- unlist(results[(s+1): (2*s)])
        ubs[i, ] <- unlist(results[(2*s +1): length(results)])
      }
      pb$tick()
    }

    opar <- graphics::par()
    rows <- floor(sqrt(length(theta) * (16/9)))
    cols <- ceiling(length(theta) / rows)
    graphics::par(mfrow = c(rows, cols), family = "serif")

    hist(lbs[, 1], main = "Confidence interval distributions for lambda hat",
         xlim = c(theta[1] - 0.4, theta[1] + 0.4), col = "green")
    hist(ubs[, 1], col = "red", add = TRUE)
    graphics::abline(v=theta[1])

    for(i in 2:length(theta)){
      hist(lbs[, i], main = paste0("Confidence interval distributions for beta", i-2, " hat"),
           xlim = c(theta[i] - 0.4, theta[i] + 0.4), col = "green")
      hist(ubs[, i], col = "red", add = TRUE)
      graphics::abline(v=theta[i])
    }
    suppressWarnings(graphics::par(opar))
  }

  return_vector <- vector()
  for(i in 1:length(theta)){
    return_vector <- append(return_vector, sum(TF_table[,i] / monte_carlo_iterations))
  }
  return(return_vector)
}

#_________________________________PRINTING______________________________________
print_as_kable <- function(eval, latex = F){
  if(latex){
    kable(eval, format = "latex" ,linesep = "")
  }else{
    kable(eval)
  }
}

#' Avaliação do modelo de regressao quantílico Chen
#' @description
#' Considerando o modelo de regressão quantílico proposto por Souza(2021) esta função
#' busca um modo prático de análisar as propriedades dos estimadores de máxima
#' verossimilhança e transformá-las em Latex.
#' @param n Inteiro que representa o número de valores aleatorios a serem gerados
#' no cálculo das propriedades do estimador.
#' @param theta Vetor de numéricos com os parâmetros. O primeiro é considerado lambda,
#' todos os outros são considerados parâmetros do modelo de regressão quantílico Chen.
#'
#' @return
#' Uma tabela formatada em Latex utilizando o pacote kable. Nela constamas seguintes propriedades do estimador:
#' Média, Viés; Erro padrão; Erro quadrado médio; Limite inferior do intervalo de
#' confiança; Limite superior do intervalo de confiança; Taxa de cobertura.
#' @export
#' @importFrom kableExtra kable
#' @import progress
#'
#' @examples
#' theta1 <- c(0.7, 2, -1, 0.8)
#' make_table(30, theta1)
#' make_table(100, theta1)
#' make_table(200, theta1)
#'
#' theta2 <- c(0.7, 1.5, -0.8, 0.5, 2, -1)
#' make_table(30, theta2)
#' make_table(100, theta2)
#' make_table(200, theta2)


make_table <- function(n, theta){
  estimacoes <- t(eval_estim(n, 5000, theta))
  ci <- ci(n, theta, 0.05)
  taxa_de_cobertura <- c(eval_ci(n, 5000, theta, 0.05))

  cenario <- cbind(estimacoes, ci, taxa_de_cobertura)

  colnames(cenario) <- c("Média", "Viés", "Erro padrão", "EQM", "LI", "LS", "TC")
  cenario = cenario[2:nrow(cenario),]
  return(print_as_kable(cenario ,latex = T))
}

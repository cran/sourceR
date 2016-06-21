## ----run_doMCMC_sim_pois, results='hide', message=FALSE, eval = FALSE, cache = T----
#  require(sourceR)
#  set.seed(63164)
#  data(sim_SA)
#  data(sim_SA_true)
#  
#  # Set priors
#  priors <- list(a = 1, r = 1, theta = c(0.01, 0.00001))
#  
#  # Run model
#  res_sim <- saBayes(formula = Human ~ Source1 + Source2 + Source3 + Source4 + Source5,
#                 time = ~Time, location = ~Location, type = ~Type,
#                 data = sim_SA$data, priors = priors,
#                 alpha_conc = 1, prev = sim_SA$prev,
#                 likelihood_dist = "pois", n_iter = 1010,
#                 mcmc_params = list(burn_in = 20, thin = 1))

## ----trace_acf_sim_data_code, dev='tikz',results='hide', eval=FALSE, cache = T----
#  ## Plot the marginal posterior for the source effect 2, at time 1, location A
#  plot(res_sim$posterior$a$time1$locationA[,"Source3"], type="l")
#  ## Plot the marginal posterior for the type effect 21
#  plot(res_sim$posterior$q[,"type21"], type="l")
#  ## Plot the marginal posterior for the relative prevalence of source effect 5,
#  ## type 17, at time 2
#  plot(res_sim$posterior$r$time2["type17","Source5",], type="l")

## ----summary_sim_pois, results='hide', message=FALSE, eval = FALSE, cache = T----
#  summary(res_sim, alpha = 0.05, thin = 1, burn_in = 0)

## ----subset_sim_pois, results='hide', message=FALSE, eval = FALSE, cache = T----
#  subset_posterior(res_sim, params = c("a", "li", "q"),
#                   t = "1", l = "B", j = c("Source2", "Source1"),
#                   i = c("47", "10"), iters = c(3:10))

## ----flatten_sim_pois, results='hide', message=FALSE, eval = FALSE, cache = T----
#  flatten(res_sim)

## ----run_doMCMC_real, results='hide', eval = FALSE, message=FALSE, cache = T----
#  data(campy)
#  set.seed(59623)
#  # remove rows with no source cases as there is no information for
#  # source attribution of human cases for these sources
#  zero_rows <- which(apply(campy[,c(2 : 7)], 1, sum) == 0)
#  campy <- campy[-zero_rows,]
#  
#  # Set priors
#  priors <- list(a = 1, r = 1, theta = c(0.01, 0.00001))
#  
#  # set prevalences
#  # Number of samples  tested  for c. jejuni, for each source.
#  tot_samples<-c(239, 196, 127, 595, 552, 192 + 332)
#  # Number of samples positive for c. jejuni, for each source.
#  pos_samples<-c(181, 113, 109, 97, 165, 24 + 62)
#  prevs <- data.frame(value = pos_samples / tot_samples,
#                      source_id = colnames(campy[, 2:7]))
#  
#  # Run model
#  # the model assumes one time and location if none are specified in saBayes
#  res_real <- saBayes(formula = Human ~ ChickenA + ChickenB + ChickenC +
#                   Bovine + Ovine + Environment,
#                 type = ~Type, data = campy, priors = priors, alpha_conc = 1,
#                 prev = prevs, likelihood_dist = "pois", n_iter = 1020,
#                 mcmc_params = list(burn_in = 20, thin = 1))


Guttmans_reliability <- function(Data) {
  ftt <- alpha(Data)
  ftt_fit <-c(ftt[["item.stats"]][["mean"]],
                  ftt[["item.stats"]][["sd"]],
                  ftt[["item.stats"]][["r.cor"]],
                  ftt[["alpha.drop"]][["std.alpha"]],
                  ftt[["alpha.drop"]][["G6(smc)"]])
  ftt_mm <- c(ftt[["total"]][["std.alpha"]],
              ftt[["total"]][["G6(smc)"]],
              ftt[["total"]][["average_r"]],
              ftt[["total"]][["mean"]],
              ftt[["total"]][["sd"]])
  
  rownames(ftt_mm) <- c("Alpha","Guttmans Lamda reliability",
                  "Average interitem correlation","Mean score","SD")
  rownames(ftt_fit) <- c("Mean ","SD ","correlation","Alpha ","Reliability")
}

#parallel factoring
factoring <- function(Data, exclude) {
  PCA_items <-rbind("degree of freedom",
                "Chi-sq","Chi-sq/df",
                "Harmonic sample size",
                "Root Mean Square",
                "Probability of the empirical chi-sq",
                "Adjusted Root Mean Square",
                "Empirical BIC","Sample size adjusted BIC",
                "fit (SSresidual vs SSoriginal values)",
                "fit applied to off diagonal elements",
                "SD of the residuals",
                "Number of factors extracted",
                "Number of observations",
                "Value of the minimised function",
                "chi-sq based on the objective function",
                "p-value of observing the chi-sq",
                "chi-sq based on the objective function/df",
                "Null model", "df for null model","chi-sq for null model",
                "chi-sq for null model/df",
                "Tucker Lewis Index of factoring reliability",
                "RMSE Approximation", "RMSE Approximation-lower",
                "RMSE Approximation-upper",
                "RMSE Approximation-confidence interval",
                "RMSE Approximation-BIC",
                "RMSE Approximation-empirical BIC",
                "Mean item complexity",
                "Kaiser Meyer Olkin Measure of Sampling Adequacy",
                "Bartlett Chi", "Barlett p-value",
                "Barlett df", "Barlett Chi/df")
  
  lll <- fa.parallel(Data, fm = 'minres', se.bars = TRUE, fa = 'fa')
  
  ll <- fa(Data, nfactors = lll[["nfact"]], rotate = "varimax", fm="minres")
  
  km <- KMO(Data)
  
  bc <- cortest.bartlett(Data, n = nrow(Data)- exclude, diag = TRUE)
  
  Values <- rbind(ll[["dof"]], ll[["chi"]], ll[["chi"]]/ll[["dof"]], 
        ll[["nh"]], ll[["rms"]], ll[["EPVAL"]], ll[["crms"]], 
        ll[["EBIC"]], ll[["ESABIC"]], ll[["fit"]], ll[["fit.off"]],
        ll[["sd"]], ll[["factors"]], ll[["n.obs"]], 
        ll[["objective"]], ll[["STATISTIC"]], ll[["PVAL"]],
        ll[["STATISTIC"]]/ll[["null.dof"]], ll[["null.model"]], 
        ll[["null.dof"]], ll[["null.chisq"]], 
        ll[["null.chisq"]]/ll[["null.dof"]], ll[["TLI"]],
        ll[["RMSEA"]][["RMSEA"]], 
        ll[["RMSEA"]][["lower"]], 
        ll[["RMSEA"]][["upper"]],
        ll[["RMSEA"]][["confidence"]],
        ll[["BIC"]], ll[["SABIC"]], mean(ll[["complexity"]]),
        km[["MSA"]], bc[["chisq"]], bc[["p.value"]], bc[["df"]],
        bc[["chisq"]]/bc[["df"]])
  rownames(Values) <- PCA_items
  Values <- data.frame(Values)
  
  #colnames(llj) <- c("Extracted PCA Parameters", "Value estimated")
  
  results <- list(#Plottted_Model = Diagramme,
                  Extrated_factors = ll,
                  Parameters_of_Extracted_factos = Values)
  return(results)
}

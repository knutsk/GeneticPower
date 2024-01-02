#' Power for Quantitative Trait with Numeric Genotype Count
#'
#' This function returns the power for a linear regression of genetic association for a quantitative trait and a numeric treatment of genotype count.
#' @param N total samples in the analysis
#' @param delta Effect size difference between homozygote risk/disease allele vs. homozygote reference allele, aka: mean(bb) - mean (AA), where 'b' is the disease allele, 'A' is the reference allele
#' @param freq allele frequency of disease (alternate) allele 'b'
#' @param sigma standard deviation of the response phenotype
#' @param NCovariates the number of additional parameters in the model (will reduce your overall degrees of freedom)
#' @param alpha desired significance level of association. Default is 0.05
#' @param numtests the number of tests to be corrected by Bonferroni adjustment before achieving 'alpha'
#' @param minh mode of inheritance:  "recessive", "additive", "dominant" (this is the same as moi=0,0.5, and 1.0, respectively). Default is "additive" if no moi is specified.
#' @param moi mode of inheritance: 0 for recessive, 0.5 for additive, 1.0 for dominant, or anywhere in between. This parameter OVER-RIDES minh, and is useful for modeling in-between moi's.
#' @param rsquared fraction of total sum-of-squares explained by model fit - OVER-RIDES delta AND sigma. Useful if the only benchmark summary statistics available is an R2 value.
#' @import stats4
#' @export
#' @examples
#' GeneticPower.QuantitativeTrait.Numeric()

GeneticPower.QuantitativeTrait.Numeric <- function(N=1000,delta=1,freq=0.15,minh='additive',
                                                   sigma=1,NCovariates=0,alpha=0.05, numtests=1,moi=NULL,rsquared=NULL) {

  alphy <- 1-(1-alpha)^(1/numtests);


  if (is.null(moi)) {
    minh <- match.arg(minh) # can use abbreviated names
    moi <- switch(minh,
                  additive   = 0.5,
                  dominant   = 1.0,
                  recessive  = 0  ) # define the mode of inheritance
  }


  N1 <- N*(1-freq)^2;
  N2 <- 2*N*freq*(1-freq);
  N3 <- N*freq^2;


  mu <- c(0,moi,1)*delta;

  #NewRsquared <- ((N1+N3)*N2*beta^2 + (N1+N2)*N3*delta^2 - 2*N2*N3*beta*delta)/(((N1+N2+N3)*sigma)^2 + (N1+N3)*N2*beta^2 + (N1+N2)*N3*delta^2 - 2*N2*N3*beta*delta);

  if (is.null(rsquared)) {

    xbar <- (N2 + N3*2)/N;
    ybar <- (N2*mu[2] + N3*mu[3])/N;

    slope <- (N1*xbar*ybar + (1-xbar)*N2*(mu[2]-ybar) + (2-xbar)*N3*(mu[3]-ybar))/(N1*xbar^2 + N2*(1-xbar)^2 + N3*(2-xbar)^2);
    a <- ybar - slope*xbar;

    fits <- a + slope*c(0,1,2);

    eNs <- c(N1,N2,N3);

    Rsqtop <- sum(eNs*(fits-ybar)^2);

    #Rsqbottom <- Rsqtop + sum((eNs-1)*sigma^2 + eNs*(mu - fits) + eNs*(mu-fits)^2);
    Rsqbottom <- Rsqtop + sum((eNs-1)*sigma^2 + 2*eNs*(mu - fits)*(fits-ybar) + eNs*(mu-fits)^2);

    NewRsquared <- Rsqtop/Rsqbottom;

  } else {

    NewRsquared <- rsquared;
  }

  power <- pf(qf(1-alphy,df1=1,df2=(N-2-NCovariates)),ncp=(NewRsquared*(N-2)/((1-NewRsquared))),df1=1,df2=(N-2-NCovariates),lower.tail=F);

  return(power);
}


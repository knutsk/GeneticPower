

#' Effect Size Required for Powering a Genetic Association for a Quantiative Trait with Numeric Genotype Count
#'
#'This function returns the required effect size (genetic effect difference between homozygote risk vs homozygote ref) to power a linear regression of genetic association for a quantitative trait and a numeric treatment of genotype count.
#' @param power the required statistical power as a decimal (e.g., for 80% power, use 0.8)
#' @param N total samples in the analysis
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
#' GeneticEffect.QuantitativeTrait.Numeric()

GeneticEffect.QuantitativeTrait.Numeric<- function(power=0.8,N=5000,freq=0.15,sigma=1,NCovariates=0,alpha=0.05,numtests=1,minh="additive",moi=0.5,rsquared=NULL) {

  A <- 1e6*sigma;
  flag <- 0;
  B <- 0;
  while(A-B > .001*sigma){
    ga <- (A+B)/2
    result <- GeneticPower.QuantitativeTrait.Numeric(N=N,delta=ga,freq=freq,minh=minh,sigma=sigma,NCovariates=NCovariates,alpha=alpha, numtests=numtests,moi=moi,rsquared=rsquared)

    po <-result;

    if (po > power) {A <- ga; flag <- 1;} else {B <- ga;}

  }


  if (flag == 1) {
    return(ga);
  } else {
    return(NA);
  }


}

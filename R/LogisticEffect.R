#' Effect Size Required for Powering a Case-Control Design with Numeric Genotype Count
#'
#' This function returns the required effect size (homozygote risk vs homozygote ref) to power a logistic regression of genetic association for a case-control/binary trait and a numeric treatment of genotype count.
#' @param power the required statistical power as a decimal (e.g., for 80% power, use 0.8)
#' @param N total samples in the anlysis, including both cases and control
#' @param fc "fraction of cases" = fraction of collected samples considered to be a disease "case"
#' @param freq allele frequency of disease (alternate) allele 'b'
#' @param kp fraction of overall HWE population that has the disease/ NOTE: This is NOT necessarily the same as the fraction of collected sample with the disease (this is the population-level prevalence)
#' @param alpha desired significance level of association. Default is 0.05
#' @param misclass.control the fraction of subjects labeled as controls who are mis-classified as controls (and are hence really cases)
#' @param misclass.case the fraction of subjects labeled as cases who are mis-classified as cases (and are hence really controls)
#' @param minh mode of inheritance:  "recessive", "additive", "dominant" (this is the same as moi=0,0.5, and 1.0, respectively). Default is "additive" if no moi is specified.
#' @param moi mode of inheritance: 0 for recessive, 0.5 for additive, 1.0 for dominant, or anywhere in between. This parameter OVER-RIDES minh, and is useful for modeling in-between moi's.
#' @param numtests the number of tests to be corrected by Bonferroni adjustment before achieving 'alpha'
#' @import stats4
#' @export
#' @examples
#' GeneticEffect.CaseControlTrait.NumericGenotype()

GeneticEffect.CaseControlTrait.NumericGenotype <- function(power,N,fc=0.5,freq,kp,alpha=5.0e-08, misclass.control = 0, misclass.case=0,minh='additive',numtests = 1, moi = NULL) {

  A <- 1e6;
  flag <- 0;
  B <- 1.0001;
  while(A/B > 1.001){
    ga <- sqrt(A*B);
    result <- GeneticPower.CaseControl.Logistic.Numeric(N=N,fc=fc,gamma=ga,freq=freq,kp=kp,alpha=alpha);
    GeneticPower.CaseControlTrait.NumericGenotype(N=N,fc=fc,gamma=ga,freq=freq,kp=kp,alpha=alpha,misclass.control=misclass.control,misclass.case=misclass.case,minh=minh,numtests = numtests, moi = moi)
    po <-result$power;

    if (po > power) {A <- ga; flag <- 1;} else {B <- ga;}

  }


  if (flag == 1) {
    return(ga);
  } else {
    return(NA);
  }


}

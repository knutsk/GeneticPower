
#' Sample Size Required for Powering a Case-Control Design with Numeric Genotype Count
#'
#' This function returns the required sample size (Total N, aka cases + controls) to power a logistic regression of genetic association for a case-control/binary trait and a numeric treatment of genotype count. You can calculate the number of cases by multiplying the total N by fc.
#' @param power the required statistical power as a decimal (e.g., for 80% power, use 0.8)
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
#' N.CaseControlTrait.NumericGenotype()

N.CaseControlTrait.NumericGenotype <- function(power=0.8,gamma=1.5,freq=0.15, kp=.1, alpha=5.0e-08, fc=0.5 ,misclass.case=0,misclass.control=0,moi=0.5, minh = "additive", numtests = 1) {

  A <- 1e9;
  flag <- 0;
  B <- 0;
  while(A-B > 5){
    ga <- (A+B)/2;

    result <- GeneticPower.CaseControlTrait.NumericGenotype(N=ga,fc=fc,gamma=gamma,freq=freq,kp=kp,alpha=alpha,misclass.control=misclass.control,misclass.case=misclass.case,minh=minh,numtests = numtests, moi = moi)

    po <-result$power;

    if (po > power) {A <- ga; flag <- 1;} else {B <- ga;}

  }


  if (flag == 1) {
    return(ga);
  } else {
    return(NA);
  }


}


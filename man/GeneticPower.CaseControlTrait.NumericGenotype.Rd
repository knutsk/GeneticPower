% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LogisticPower.R
\name{GeneticPower.CaseControlTrait.NumericGenotype}
\alias{GeneticPower.CaseControlTrait.NumericGenotype}
\title{Power for Case-Control Design with Numeric Genotype Count}
\usage{
GeneticPower.CaseControlTrait.NumericGenotype(
  N,
  gamma,
  freq,
  kp,
  alpha = 0.05,
  fc = 0.5,
  misclass.control = 0,
  misclass.case = 0,
  minh = "additive",
  numtests = 1,
  moi = NULL
)
}
\arguments{
\item{N}{total samples in the anlysis, including both cases and control}

\item{gamma}{odds ratio between disease-allele homozygote (bb) and reference allele homozygote (AA), i.e. odds(D|bb)/odds(D|AA).}

\item{freq}{allele frequency of disease (alternate) allele 'b'}

\item{kp}{fraction of overall HWE population that has the disease/ NOTE: This is NOT necessarily the same as the fraction of collected sample with the disease (this is the population-level prevalence)}

\item{alpha}{desired significance level of association. Default is 0.05}

\item{fc}{"fraction of cases" = fraction of collected samples considered to be a disease "case"}

\item{misclass.control}{the fraction of subjects labeled as controls who are mis-classified as controls (and are hence really cases)}

\item{misclass.case}{the fraction of subjects labeled as cases who are mis-classified as cases (and are hence really controls)}

\item{minh}{mode of inheritance:  "recessive", "additive", "dominant" (this is the same as moi=0,0.5, and 1.0, respectively). Default is "additive" if no moi is specified.}

\item{numtests}{the number of tests to be corrected by Bonferroni adjustment before achieving 'alpha'}

\item{moi}{mode of inheritance: 0 for recessive, 0.5 for additive, 1.0 for dominant, or anywhere in between. This parameter OVER-RIDES minh, and is useful for modeling in-between moi's.}
}
\description{
This function returns the power for a logistic regression of genetic association for a case-control/binary trait and a numeric treatment of genotype count.
}
\examples{
GeneticPower.CaseControlTrait.NumericGenotype()
}

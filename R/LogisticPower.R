#' Power for Case-Control Design with Numeric Genotype Count
#'
#' This function returns the power for a logistic regression of genetic association for a case-control/binary trait and a numeric treatment of genotype count.
#' @param N total samples in the anlysis, including both cases and control
#' @param gamma odds ratio between disease-allele homozygote (bb) and reference allele homozygote (AA), i.e. odds(D|bb)/odds(D|AA).
#' @param freq allele frequency of disease (alternate) allele 'b'
#' @param kp fraction of overall HWE population that has the disease/ NOTE: This is NOT necessarily the same as the fraction of collected sample with the disease (this is the population-level prevalence)
#' @param alpha desired significance level of association. Default is 0.05
#' @param fc "fraction of cases" = fraction of collected samples considered to be a disease "case"
#' @param misclass.control the fraction of subjects labeled as controls who are mis-classified as controls (and are hence really cases)
#' @param misclass.case the fraction of subjects labeled as cases who are mis-classified as cases (and are hence really controls)
#' @param minh mode of inheritance:  "recessive", "additive", "dominant" (this is the same as moi=0,0.5, and 1.0, respectively). Default is "additive" if no moi is specified.
#' @param moi mode of inheritance: 0 for recessive, 0.5 for additive, 1.0 for dominant, or anywhere in between. This parameter OVER-RIDES minh, and is useful for modeling in-between moi's.
#' @param numtests the number of tests to be corrected by Bonferroni adjustment before achieving 'alpha'
#' @import stats4
#' @export
#' @examples
#' GeneticPower.CaseControlTrait.NumericGenotype()

GeneticPower.CaseControlTrait.NumericGenotype <- function(N, gamma, freq, kp, alpha=.05, fc=0.5 ,misclass.control=0,misclass.case=0,minh='additive',numtests = 1, moi = NULL)
{

  if (is.null(moi)) {
    minh <- match.arg(minh) # can use abbreviated names
    moi <- switch(minh,
                  additive   = 0.5,
                  dominant   = 1.0,
                  recessive  = 0  ) # define the mode of inheritance
  }


  if ((moi < 0) || (moi > 1)) {
    print("invalid inheritance mode: must be between 0 and 1");
    return(NA);
  }

  minh <- moi;

  alphy <- 1-(1-alpha)^(1/numtests);
  pCntrl <- 1 - fc;

  p <- freq;
  q <- 1 - p # frequency of 'A' allele
  pd <- kp;  #disease prevalence in population

  #redundant names allow easy re-use of code...

  fhw <- c(p^2, 2*p*q, q^2) # frequencies of genotype according to HWE rule

  #####here, we calculate the pr(D|AA) based on kp and odds ratio

  gp <- (gamma-1)*((gamma^minh) - 1);      #gp and gs these mean nothing, they just appear a lot in the calculations...
  gs <- (gamma-1) + ((gamma^minh) - 1);

  coff.x3 <- gp*(1-freq)^2;
  coff.x2 <- (1-freq)*(1-freq)*gs + 2*freq*(1-freq)*(gamma^minh)*(gamma-1) + freq*freq*gamma*((gamma^minh)-1) - gp*pd;
  coff.x <- (1-freq)*(1-freq) + (gamma^minh)*2*freq*(1-freq) + gamma*freq*freq - gs*pd;
  coff.1 <- pd*(-1);

  ##coefficients of cubic in 'x' where x= Pr(Disease|AA), A=non-disease allele.

  coff.vect <- c(coff.1,coff.x,coff.x2,coff.x3);
  roots <- polyroot(coff.vect);

  irealposprob <- (abs(roots - abs(roots)) < 1.0e-10) & (abs(roots) <= 1);
  useroot <- abs(roots[irealposprob]);
  pii <- useroot;   #pii = pr(Disease | AA)
  rootcheck <- pd - freq*freq*(gamma*pii/(1+(gamma-1)*pii)) - 2*freq*(1-freq)*((gamma^minh)*pii/(1+(gamma^minh-1)*pii)) - pii*(1-freq)^2;

  f.mod <- c(gamma/(1+(gamma-1)*pii),(gamma^minh)/(1+(gamma^minh-1)*pii),1);

  ####Do NOt panic!!! f.mod is defined in reverse order, so yes, index 2 is for the heterozyogte

  pG0givenSick <- pii*(1-freq)*(1-freq)/pd;
  pG1givenSick <- pii*f.mod[2]*2*freq*(1-freq)/pd;
  pG2givenSick <- pii*f.mod[1]*freq*freq/pd;

  pG0givenWell <- (1-pii)*(1-freq)*(1-freq)/(1-pd);
  pG1givenWell <- (1-pii*f.mod[2])*2*freq*(1-freq)/(1-pd);
  pG2givenWell <- (1-pii*f.mod[1])*freq*freq/(1-pd);

  expSickG0 <- N*fc*(1-misclass.case)*pG0givenSick;
  expNormG0 <- N*pCntrl*(1-misclass.control)*pG0givenWell;
  expSickMarkedNormG0 <- N*pCntrl*(misclass.control)*pG0givenSick;
  expNormMarkedSickG0 <- N*fc*(misclass.case)*pG0givenWell;

  expSickG1 <- N*fc*(1-misclass.case)*pG1givenSick;
  expNormG1 <- N*pCntrl*(1-misclass.control)*pG1givenWell;
  expSickMarkedNormG1 <- N*pCntrl*(misclass.control)*pG1givenSick;
  expNormMarkedSickG1 <- N*fc*(misclass.case)*pG1givenWell;

  expSickG2 <- N*fc*(1-misclass.case)*pG2givenSick;
  expNormG2 <- N*pCntrl*(1-misclass.control)*pG2givenWell;
  expSickMarkedNormG2 <- N*pCntrl*(misclass.control)*pG2givenSick;
  expNormMarkedSickG2 <- N*fc*(misclass.case)*pG2givenWell;

  ####get parameter fits for maximum liklihood

  optim <- mle(LogisticLogLikeHood,start=list(alpha=fc,beta=0),fixed=list(Nd0=expSickG0+expNormMarkedSickG0,Nd1=expSickG1+expNormMarkedSickG1,Nd2=expSickG2+expNormMarkedSickG2,Nh0=expNormG0+expSickMarkedNormG0,Nh1=expNormG1+expSickMarkedNormG1,Nh2=expNormG2+expSickMarkedNormG2));
  fitted.intercept <- attributes(optim)$coef['alpha'];
  fitted.slope <- attributes(optim)$coef['beta'];

  ####create large idealized dataset to get best maximum likellihood estimates of linear fit for non-centrality parameter
  #
  #    Inflate <- round(100000/N); ##should get a pretty good estimate from 100,000 samples, no?
  #
  #    expResponse <- c(rep(1,round(Inflate*expSickG0)+round(Inflate*expSickG1)+round(Inflate*expSickG2)),rep(0,round(Inflate*expNormG0)+round(Inflate*expNormG1)+round(Inflate*expNormG2)));
  #    expGtypes <- c(rep(0,round(Inflate*expSickG0)),rep(1.0,round(Inflate*expSickG1)),rep(2.0,round(Inflate*expSickG2)),rep(0,round(Inflate*expNormG0)),rep(1.0,round(Inflate*expNormG1)),rep(2.0,round(Inflate*expNormG2)));
  #
  #    IdealGLM <- glm(as.factor(expResponse) ~ expGtypes,family=binomial);
  #
  #    fitted.intercept <- as.numeric(IdealGLM$coefficients[1]);
  #    fitted.slope <- as.numeric(IdealGLM$coefficients[2]);

  npDiseaseg0 <- exp(fitted.intercept)/(1+exp(fitted.intercept));
  npDiseaseg1 <- exp(fitted.intercept+fitted.slope)/(1+exp(fitted.intercept+fitted.slope));
  npDiseaseg2 <- exp(fitted.intercept+2*fitted.slope)/(1+exp(fitted.intercept+2*fitted.slope));

  ####set up the Wald Statistic ...

  #     expNormG0 <- N*pCntrl*pG0givenWell;
  #     expNormG1 <- N*pCntrl*pG1givenWell;
  #     expNormG2 <- N*pCntrl*pG2givenWell;

  expNormG0 <- expNormG0 + expSickMarkedNormG0;
  expNormG1 <- expNormG1 + expSickMarkedNormG1;
  expNormG2 <- expNormG2 + expSickMarkedNormG2;

  expSickG0 <- expSickG0 + expNormMarkedSickG0;
  expSickG1 <- expSickG1 + expNormMarkedSickG1;
  expSickG2 <- expSickG2 + expNormMarkedSickG2;

  LogNullLHood <- N*fc*log(fc) + N*pCntrl*log(pCntrl);
  LogModelLHood <- expSickG0*log(npDiseaseg0)+expNormG0*log(1-npDiseaseg0) + expSickG1*log(npDiseaseg1)+expNormG1*log(1-npDiseaseg1) + expSickG2*log(npDiseaseg2)+expNormG2*log(1-npDiseaseg2);

  lambda <- -2*(LogNullLHood-LogModelLHood);

  return(data.frame(cbind(power=pchisq(qchisq(1-alphy, df=1), df=1, ncp=lambda, lower.tail=F),pii=pii,rootcheck=rootcheck,RR.althomo=f.mod[1],RR.het=f.mod[2])));

}


#' Log Likelihood Function for Logistic Model
#'
#' This function returns the negative log likelihood for logistic regression
#' @param Nd0 Number of disease subjects with 0 risk alleles
#' @param Nd1 Number of disease subjects with 1 risk allele
#' @param Nd2 Number of disease subjects with 2 risk alleles
#' @param Nh0 Number of healthy subjects with 0 risk alleles
#' @param Nh1 Number of healthy subjects with 1 risk allele
#' @param Nh2 Number of healthy subjects with 2 risk alleles
#' @param alpha Intercept
#' @param beta Odds Ratio for Genotype Effect
#' @export
#' @examples
#' LogisticLogLikeHood()

LogisticLogLikeHood <- function(Nd0,Nd1,Nd2,Nh0,Nh1,Nh2,alpha=0.5,beta=0.01) {
  res <- Nd0*alpha + Nd1*(alpha+beta) + Nd2*(alpha + 2*beta) - (Nd0+Nh0)*log(1+exp(alpha)) - (Nd1+Nh1)*log(1+exp(alpha+beta)) - (Nd2+Nh2)*log(1+exp(alpha+2*beta));
  return(-res);
}

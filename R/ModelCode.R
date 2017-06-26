#' nmadasmodel
#' @description  Specify the copula based bivariate beta-binomial or alternatively logistic-binomial distribution to fit to the diagnostic data.
#' @param marginals Use normal marginals on the logit transformed sensitivity and specificity or the beta marginals.
#' When marginals = 'normal' the following model is fitted:
#' \deqn{Y_{ijk} ~ bin(\pi_{ijk}, N_{ijk})}
#' \deqn{logit(\pi_{ijk}) = \mu_{jk} + \eta_{ij} + \delta_{ijk}}
#' \deqn{(\eta_{i1}, \eta_{i2})' ~ N_2(0, \Sigma)}
#' \deqn{\Sigma[j,j] = \sigma[j]^2, \Sigma[1,2] = \Sigma[2,1] = \rho*\sigma[1]*\sigma[2]}
#' \deqn{\delta_{ijk} ~ N(0, \tau_{jk})}
#' @param fullcov Logical for full (TRUE) or reduced (default) variance-covariance matrix. The reduction simplifies
#' the variance-covariance matrix by specifying that
#' \deqn{\delta_{ijk} ~ N(0, \tau_{j})}.
#' @param prior.lmu A text specifying the prior distribution for \eqn{\mu}  parameters. The default is \code{"normal(0, 5)"}.
#' @param prior.tau A text specifying the prior distribution for \eqn{\tau}  parameters. The default is \code{"cauchy(0, 2.5)"}.
#' @param prior.sigma A text specifying the prior distribution for \eqn{\sigma}  parameters. The default is \code{"cauchy(0, 2.5)"}.
#' @param prior.rho A text specifying the prior distribution for \eqn{\rho}  parameter. The default is \code{"lkj(2.0)"}.
#' For more details on specifications of prior distributions see \href{http://mc-stan.org/documentation/}{Stan documentation}.
#' When marginals = 'beta' the following model is fitted:
#' \deqn{Y_{ijk} ~ bin(\pi_{ijk}, N_{ijk})}
#' \deqn{\pi_{i1k}, \pi_{i2k} ~ f(\pi_{i1k})*f(\pi_{i2k})*copula(F(\pi_{i1k}), F(\pi_{i2k}), \omega_k)}
#' where f and F are the probability density and cumulative distribution function of a beta distribution with parameters
#' \deqn{\alpha_{jk}} and \deqn{\beta_{jk}} specified as follows:
#' \deqn{\alpha_{jk} = \mu_{jk}*\frac{1- \theta_j*\delta_{jk}}{\theta_j*\delta_{jk}}}
#' \deqn{\beta_{jk} = (1 - \mu_{jk})*\frac{1- \theta_j*\delta_{jk}}{\theta_j*\delta_{jk}}}
#' Here \deqn{\mu_{jk}} is the mean sensitivity \code{j = 1} and specificity \code{j = 2},
#'  \deqn{\omega_k} captures the correlation between sensitivity and specificity in test k,
#'  \deqn{\theta_j} captures the common overdispersion among the sensitivities \code{j = 1} and specificities \code{j = 2},
#'  and \deqn{\delta_jk} captures the test specific extra variability.
#'  The hyper parameters \eqn{\mu}, \eqn{\theta} and \eqn{\delta} are given beta/uniform priors since they are in the (0,1) interval.
#'  The prior distribution of \eqn{\omega} depends on the copula.
#' @param copula Name of the copula function used to model the correlation between sensitivity and specificty. This
#' requires that marginals = 'beta' be specificied.
#' This is a string naming the copula function. The choices are "fgm", "frank", "gauss", "c90" and "c270".
#' @param p.omega The prior distribution of the \eqn{\omega} parameters. This prior distribution depend on the
#'  specified copula. The defualt is \code{"uniform(-1, 1)"} for the gaussian and fgm copula,
#'  \code{"normal(0, 5)"} for the frank copula, ad \code{"cauchy(0, 2.5)"} for the c90 and c270 copula.
#' @examples
#' model1 <-  nmamodel()
#'
#' model2 <- nmamodel(copula = 'fgm')
#'
#' model3 <-  nmamodel(marginals = 'normal')
#'
#' @references {Agresti A (2002). Categorical Data Analysis. John Wiley & Sons, Inc.}
#' @references {Clayton DG (1978). A model for Association in Bivariate Life Tables and its Application in
#' Epidemiological Studies of Familial Tendency in Chronic Disease Incidence. Biometrika,65(1), 141-151.}
#' @references {Frank MJ (1979). On The Simultaneous Associativity of F(x, y) and x + y - F(x, y). Aequationes Mathematicae, pp. 194-226.}
#' @references {Farlie DGJ (1960). The Performance of Some Correlation Coefficients for a General Bivariate
#' Distribution. Biometrika, 47, 307-323.}
#' @references {Gumbel EJ (1960). Bivariate Exponential Distributions. Journal of the American Statistical Association, 55, 698-707.}
#' @references {Meyer C (2013). The Bivariate Normal Copula. Communications in Statistics - Theory and Methods, 42(13), 2402-2422.}
#' @references {Morgenstern D (1956). Einfache Beispiele Zweidimensionaler Verteilungen. Mitteilungsblatt furMathematische Statistik, 8, 23 - 235.}
#' @references {Sklar A (1959). Fonctions de Repartition a n Dimensions et Leurs Marges. Publications de l'Institut de Statistique de L'Universite de Paris, 8, 229-231.}
#' @export
#' @author Victoria N Nyaga
#' @return An object of nmamodel class.

nmadasmodel <- function(
				marginals = "beta",
				copula = "frank",
				p.omega = NULL,
				fullcov = FALSE,
				prior.lmu = "normal(0, 5)",
				prior.tau ="cauchy(0, 2.5)",
				prior.sigma = " cauchy(0, 2.5)",
				prior.rho = "lkj_corr(2.0)"
                  ){
  if (marginals == "beta"){
    GAUSS <- "
    functions{
    	real gaussian_lpdf(matrix p_i, vector alpha, vector beta, real omega){
    		real f1;
    		real f2;
    		vector[rows(p_i)] f3;
    		vector[rows(p_i)] f4;
    		int r;

    		r = rows(p_i);
    		f1 = beta_lpdf(col(p_i, 1)| alpha[1], beta[1]);
    		f2 = beta_lpdf(col(p_i, 2)| alpha[2], beta[2]);
    		for (i in 1:r){
    			f3[i] = (1/sqrt(1 - omega^2))*exp((2*omega*inv_Phi(beta_cdf(p_i[i, 1], alpha[1], beta[1]))*inv_Phi(beta_cdf(p_i[i, 2], alpha[2], beta[2])) -
    			omega^2*(inv_Phi(beta_cdf(p_i[i, 1], alpha[1], beta[1]))^2 + inv_Phi(beta_cdf(p_i[i, 2], alpha[2], beta[2]))^2))/
    			(2*(1 - omega^2)));
    		}
    		return (f1 + f2 + sum(log(f3)));
    	}
    }"

    C90 <- "
    functions{
    	real clayton90_lpdf(matrix p_i, vector alpha, vector beta, real omega){
    		real f1;
    		real f2;
    		vector[rows(p_i)] f3;
    		vector[rows(p_i)] f4;
    		int r;
    		real powr;

    		r = rows(p_i);
    		f1 = beta_lpdf(col(p_i, 1)| alpha[1], beta[1]);
    		f2 = beta_lpdf(col(p_i, 2)| alpha[2], beta[2]);
    		powr = -(2*omega + 1)/(omega);

    		for (i in 1:r){
    			f3[i] = (1 + omega)*((1 - beta_cdf(p_i[i,1], alpha[1], beta[1]))^(-(1 + omega)))*(beta_cdf(p_i[i,2], alpha[2], beta[2])^(-(1 + omega)))*
    			(((1 - beta_cdf(p_i[i,1], alpha[1], beta[1]))^(-omega) + beta_cdf(p_i[i,2], alpha[2], beta[2])^(-omega) - 1)^powr);
    		}
    		return (f1 + f2 + sum(log(f3)));
    	}
    }
    "

    C270 <- "
    functions{
    	real clayton270_lpdf(matrix p_i, vector alpha, vector beta, real omega){
    		real f1;
    		real f2;
    		vector[rows(p_i)] f3;
    		vector[rows(p_i)] f4;
    		int r;
    		real powr;

    		r = rows(p_i);
    		f1 = beta_lpdf(col(p_i, 1)| alpha[1], beta[1]);
    		f2 = beta_lpdf(col(p_i, 2)| alpha[2], beta[2]);
    		powr = -(2*omega + 1)/(omega);

    		for (i in 1:r){
    			f3[i] = (1 + omega)*(beta_cdf(p_i[i,1], alpha[1], beta[1])^(-(1 + omega)))*((1 - beta_cdf(p_i[i,2], alpha[2], beta[2]))^(-(1 + omega)))*
    			((beta_cdf(p_i[i,1], alpha[1], beta[1])^(-omega) + (1 - beta_cdf(p_i[i,2], alpha[2], beta[2]))^(-omega) - 1)^powr);
    		}
    		return (f1 + f2 + sum(log(f3)));
    		}
    }
    "

    FRANK <- "
    functions{
        real frank_lpdf(matrix p_i, vector alpha, vector beta, real omega){
            real f1;
            real f2;
            int r;
            vector[rows(p_i)] f3;
            vector[rows(p_i)] f4;

            r = rows(p_i);

        	f1 = beta_lpdf(col(p_i, 1)| alpha[1], beta[1]);
            f2 = beta_lpdf(col(p_i, 2)| alpha[2], beta[2]);
            for (i in 1:r){
                f3[i] = (omega*(1 - exp(-omega))*exp(-omega*(beta_cdf(p_i[i,1], alpha[1], beta[1]) + beta_cdf(p_i[i,2], alpha[2], beta[2]))));
                f4[i] = ((1 - exp(-omega)) - (1 - exp(-omega*beta_cdf(p_i[i,1], alpha[1], beta[1])))*(1 - exp(-omega*beta_cdf(p_i[i,2], alpha[2], beta[2]))));
            }
            return (f1 + f2 + sum(log((f3)./((f4).*(f4)))));
        }
    }"

    FGM <-  "
    functions{
    	real fgm_lpdf(matrix p_i, vector alpha, vector beta, real omega){
    	real f1;
    	real f2;
    	vector[rows(p_i)] f3;
        int r;

        r = rows(p_i);
    	f1 = beta_lpdf(col(p_i, 1)| alpha[1], beta[1]);
    	f2 = beta_lpdf(col(p_i, 2)| alpha[2], beta[2]);

    	for (i in 1:r){
    		f3[i] = log(1 + omega*(1 - 2*beta_cdf(p_i[i,1], alpha[1], beta[1]))*(1 - 2*beta_cdf(p_i[i,2], alpha[2], beta[2])));
    	}
    	return (f1 + f2 + sum(f3));
    	}
    }"

    if (copula=="gauss"){
    		copkies <- GAUSS
    	} else if (copula=="frank") {
    		 copkies <- FRANK
    	} else if (copula=="fgm") {
    		 copkies <- FGM
    	} else if (copula=="c90") {
    		 copkies <- C90
    	} else if (copula=="c270") {
    		 copkies <- C270
    	} else {
    		stop("Invalid copula chosen")}
  }

  d.ata <- "
  data{
      int N;
      int Nt;
      int Ns;
      int TP[N];
      int Dis[N];
      int TN[N];
      int NDis[N];
      int Study[N];
      int Test[N];
      int CIndex;
  }"
  if (marginals != "beta") {
  t.data <- "
  transformed data{
      vector[2] zero;
      int<lower=1, upper=Nt> d[2*Nt];

      zero[1] = 0;
      zero[2] = 0;

      for (r in 1:(2*Nt)){
          d[r] = (r <= Nt) ? r : (r - Nt);
      }
  }"
  }

  if (marginals == "beta") {
    beta.parameters <- "
      parameters{
          vector<lower=0, upper=1>[2] MU[Nt]; //mu_jk
          vector<lower=0, upper=1>[2] delta[Nt]; //delta_jk
          vector<lower=0, upper=1>[2] theta; //theta_j
      	  matrix<lower=0, upper=1>[Ns, 2] p_i[Nt]; //pi_ijk"
	if (copula %in% c("gauss", "fgm")){
  		omega.parameter <- "\n\t\tvector<lower = -1, upper = 1>[Nt] omega; //omega_k \n}"
  	} else if (copula == "frank"){
  		omega.parameter <- "\n\t\tvector[Nt] omega; //omega_k \n}"
  	} else if (copula %in% c("c90", "c270")){
  		omega.parameter <- "\n\t\tvector<lower = 0>[Nt] omega; //omega_k \n}"
  	}
  }
  else {
    normal.parameters <- "
      parameters{
          matrix[2, Nt] logitmu;
          matrix[Ns, 2] nu;
          matrix[Ns, 2] delta[Nt];
      	  vector<lower=0>[2] sigmab;
      	  corr_matrix[2] rhob;"
    	if(fullcov)
    	  taup <- "\tvector<lower=0>[Nt] tau[2];\n}"
    	else
    	  taup <- "\n\tvector<lower=0>[2] tau;\n}"
  }

  t.parameters = "
  transformed parameters{
      vector[2] RR[Nt];
      vector[Nt] S;
      matrix[Nt, Nt] A;
      matrix[Nt, Nt] B;
      matrix[Nt, Nt] C;
  "
  if (marginals == "beta") {
  beta.t.parameters <- "
  	vector<lower=0>[2] alpha[Nt]; //alpha_jk
  	vector<lower=0>[2] beta[Nt]; //beta_jk
  "
  }
  else{
  normal.t.parameters <- "
      vector[2] MU[Nt];
      matrix[Ns, 2] p_i[Nt];
  	vector<lower=0>[2] sigmabsq;
      matrix[2*Nt, 2*Nt] sigmawsq;"

      if(fullcov)
        tautp <-"\nvector<lower=0>[Nt] tausq[2];\n"
  		else
		    tautp <-"\n\tvector<lower=0>[2] tausq;\n"

  if (fullcov)
    tausq <- "\n\tfor (j in 1:2){
  			tausq[j] = (tau[j]).*(tau[j]);
  		}"
  	else
  	  tausq <- "\n\ttausq = (tau).*(tau);"
  }
  if (marginals == "beta") {
  beta.tp.block <- "
  	for(k in 1:Nt){
  		alpha[k] = (MU[k]).*((1 - (theta).*(delta[k]))./((theta).*(delta[k])));
  		beta[k] = (1 - MU[k]).*((1 - (theta).*(delta[k]))./((theta).*(delta[k])));
  	}
  "
  }
  else {
  normal.tp.block1 <-  "
      for (i in 1:Ns){
          for (j in 1:2){
              for (k in 1:Nt)
                  p_i[k][i,j] = inv_logit(logitmu[j,k] +  nu[i, j] + delta[k][i,j]);
          }
      }

      for (k in 1:Nt){
          for (j in 1:2){
              MU[k][j] = mean(col(p_i[k], j));
          }
      }

  	sigmabsq = (sigmab).*(sigmab);

      for (r in 1:(2*Nt)){
          for (c in 1:(2*Nt)){\n"
  sigmawsq <- if (fullcov)
  			  "\t\tif ((r <= Nt && c <= Nt) && (r == c)) sigmawsq[r,c] = sigmabsq[1] + tausq[1][d[r]] ;
  			  if  ((r > Nt && c > Nt) && (r == c)) sigmawsq[r,c] = sigmabsq[2] + tausq[2][d[r]] ;"
  		  else
  		  "\t\tif ((r <= Nt && c <= Nt) && (r == c)) sigmawsq[r,c] = sigmabsq[1] + tausq[1];
  		  if  ((r > Nt && c > Nt) && (r == c)) sigmawsq[r,c] = sigmabsq[2] + tausq[2];"

  normal.tp.block2 <- "\n\t\tif  ((r <= Nt && c <= Nt) && (r != c)) sigmawsq[r,c] = sigmabsq[1];
  		if  ((r > Nt && c > Nt) && (r != c)) sigmawsq[r,c] = sigmabsq[2];
            \tif ((r <= Nt && c > Nt) || (r > Nt && c <= Nt)) sigmawsq[r,c] = sigmab[1]*sigmab[2]*rhob[1,2];
          }
      }
  "
  }
  tp.block <- "
	for (k in 1:Nt){
		RR[k] = MU[k]./MU[CIndex];
	}
      for (l in 1:Nt){
          for(m in 1:Nt){
              A[l, m] = ((MU[l][1] > MU[m][1]) && (MU[l][2] > MU[m][2]))? 1 : 0;
              B[l, m] = ((MU[l][1] < MU[m][1]) && (MU[l][2] < MU[m][2]))? 1 : 0;
              C[l, m] = ((MU[l][1] == MU[m][1]) && (MU[l][2] == MU[m][2]))? 1 : 0;
          }

          S[l] = (2*sum(row(A, l)) + sum(row(C, l)))/(2*sum(row(B, l)) + sum(row(C, l)));
      }
  	}"
  if (marginals == "beta") {
  beta.model <-  "
  model{
  	//Priors
      for (k in 1:Nt){
          MU[k] ~ uniform(0, 1);
          delta[k] ~ uniform(0, 1);
      }
      theta ~ uniform(0, 1);"

  	if (is.null(p.omega)){
  		if (copula=="gauss"){
  			prior.omega <- "\n\t omega ~ uniform(-1, 1);\n"
  		} else if (copula=="frank"){
  			prior.omega <- "\n\t omega ~ normal(0, 5);\n"
  		} else if (copula=="fgm"){
  			prior.omega <- "\n\t omega ~ uniform(-1, 1);\n"
  		} else if (copula=="c90" | copula=="c270"){
  			prior.omega <- "\n\t omega ~ cauchy(0, 2.5);\n"
  		}
  	}
  	else{
  		prior.omega <- paste("\n\t omega ~ ", p.omega,";\n", sep="")
  	}

  	if (copula=="gauss"){
  		prior.pi <- "\n\t for (k in 1:Nt) \n\t\t p_i[k] ~ gaussian(alpha[k], beta[k], omega[k]);\n"
  	} else if (copula=="frank"){
  		prior.pi <- "\n\t for (k in 1:Nt) \n\t\t p_i[k] ~ frank(alpha[k], beta[k], omega[k]);\n"
  	} else if (copula=="fgm"){
  		prior.pi <- "\n\t for (k in 1:Nt) \n\t\t p_i[k] ~ fgm(alpha[k], beta[k], omega[k]);\n"
  	} else if (copula=="c90"){
  		prior.pi <- "\n\t for (k in 1:Nt) \n\t\t p_i[k] ~ clayton90(alpha[k], beta[k], omega[k]);\n"
  	} else if(copula=="c270"){
  		prior.pi <- "\n\t for (k in 1:Nt) \n\t\t p_i[k] ~ clayton270(alpha[k], beta[k], omega[k]);\n"
  	}
  }
  else {
  normal.model1 <- "
  model{
  	//Priors
      for (i in 1:Ns){
          nu[i] ~ multi_normal(zero, quad_form_diag(rhob, sigmab));
      }
      for (j in 1:2){\n"

  normal.priors <- paste(
  					"\t\tlogitmu[j] ~ ", prior.lmu, ";",
  					"\n\t\ttau[j] ~ ", prior.tau, ";\n\t}",
  					"\n\t\tsigmab ~ ", prior.sigma, ";\n",
                      "\t\trhob ~ ", prior.rho, ";\n",
  					sep='')
  normal.model2 <- "
      for (i in 1:Ns){
          for (j in 1:2){
              for (k in 1:Nt)\n"
  deltap <- 	if (fullcov)
  				"\t\t\tdelta[k][i,j] ~ normal(0, tau[j][k]);\n\t\t}\n\t\t}"
  			else
  			"\t\t\tdelta[k][i,j] ~ normal(0, tau[j]);\n}\n}"
  }
  last.part <-  "
      for (n in 1:N){
          TP[n] ~ binomial(Dis[n], p_i[Test[n]][Study[n], 1]);
          TN[n] ~ binomial(NDis[n], p_i[Test[n]][Study[n], 2]);
      }

  }
  generated quantities{
	vector[2*N] loglik;

      for (n in 1:N)
          loglik[n] = binomial_lpmf(TN[n]| NDis[n], p_i[Test[n]][Study[n], 1]);

      for (n in (N+1):(2*N))
          loglik[n] = binomial_lpmf(TN[n-N]| NDis[n-N], p_i[Test[n-N]][Study[n-N], 2]);

  }
  "

  if (marginals=="normal") {
  	model <- paste(
  		d.ata,
  		t.data,
  		normal.parameters,
  		taup,
  		t.parameters,
  		normal.t.parameters,
  		tautp,
  		tausq,
  		normal.tp.block1,
  		sigmawsq,
  		normal.tp.block2,
  		tp.block,
  		normal.model1,
  		normal.priors,
  		normal.model2,
  		deltap,
  		last.part,
  		sep='')
  					}
  					else{
  	model <- paste(
  		copkies,
  		d.ata,
  		beta.parameters,
		  omega.parameter,
  		t.parameters,
  		beta.t.parameters,
  		beta.tp.block,
  		tp.block,
  		beta.model,
  		prior.omega,
  		prior.pi,
  		last.part,
  		sep = "")
  					}

  out <- new("nmadasmodel",
             model = model)
  out
}

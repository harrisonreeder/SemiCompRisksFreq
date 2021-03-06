// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::vec getCommonVec(const arma::vec& delta1,const arma::vec& delta2,const arma::vec& AVec, double h){
  //computes a quantity found in many gradient components of the illness-death model
	return (delta1 + delta2 + exp(-h)) / (1 + exp(h) * AVec);
}

/*******
 Piecewise Constant
 *******/

// [[Rcpp::export]]
double nlogLikPW_uni(const arma::vec& para,
                     const arma::vec& y, const arma::vec& delta, const arma::mat& X,
                     const arma::mat& basis, const arma::mat& dbasis){

  //define constants
  int p0 = basis.n_cols;
  int p1 = X.n_cols;
  int n = X.n_rows;

  arma::vec phi = para(arma::span(0, p0-1));
  //define linear predictors
  arma::vec eta(n,arma::fill::zeros);
  if(p1 > 0){
    eta = X * para(arma::span(p0, p0+p1-1));
  }
  //define splines
  arma::vec loglambda0 = dbasis * phi;
  //under left truncation, "basis" is actually difference of
  //bases of y and yL, so Lambda0 represents correct Lambda0 - Lambda0L
  arma::vec Lambda0 = basis * arma::exp(phi);

  //delta * (log lambda) + log(survival)
  double obj_val = arma::accu( delta % (loglambda0 + eta) - Lambda0 % arma::exp(eta)  );

  return(-obj_val);
}

// [[Rcpp::export]]
arma::vec ngradPW_uni(const arma::vec& para,
                      const arma::vec& y, const arma::vec& delta, const arma::mat& X,
                      const arma::mat& basis, const arma::mat& dbasis){

  //define constants
  int p0 = basis.n_cols;
  int p1 = X.n_cols;
  int n = X.n_rows;

  arma::vec phi = para(arma::span(0, p0-1));
  //define linear predictors
  arma::vec eta(n, arma::fill::zeros);
  if(p1 > 0){
    eta = X * para(arma::span(p0, p0+p1-1));
  }
  //under left truncation, "basis" is actually difference of
  //bases of y and yL, so Lambda0 represents correct Lambda0 - Lambda0L, etc.
  arma::vec Lambda0 = basis * arma::exp(phi);
  arma::mat BLambda0 = basis.each_row() % arma::exp(phi).t();

  arma::vec temp_scorevec(p0+p1, arma::fill::zeros);
  temp_scorevec(arma::span(0, p0-1)) = dbasis.t() * delta - BLambda0.t() * arma::exp(eta);
  if(p1 > 0){
    temp_scorevec(arma::span(p0,p0+p1-1)) = X.t() * ( delta - Lambda0 % arma::exp(eta) );
  }
  return(-temp_scorevec);
}

// [[Rcpp::export]]
double nlogLikPW_ID(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
                     const arma::vec& delta1, const arma::vec& delta2,
                     const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                     const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
                     const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3,
                     const int frailty_ind){
  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
  int n = X1.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));

  double h;
  if(frailty_ind==1){
    h = para(p01+p02+p03);
  }
  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
  }

  arma::vec Lambda01 = basis1 * arma::exp(phi1);
  arma::vec Lambda02 = basis2 * arma::exp(phi2);
  //for PW model, basis and dbasis already account for Markov vs. semi-Markov
  arma::vec Lambda03 = basis3 * arma::exp(phi3);
  arma::vec loglambda01 = dbasis1 * phi1;
  arma::vec loglambda02 = dbasis2 * phi2;
  arma::vec loglambda03 = dbasis3 * phi3;

  arma::vec AVec;
  AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);

  double obj_val;
  if(frailty_ind==1){
    obj_val = arma::accu(  delta1 % ( loglambda01 + eta1 )
                             + (1-delta1) % delta2 % ( loglambda02 + eta2 )
                             + delta1 % delta2 % ( loglambda03 + eta3 + log1p( exp(h) )  )
                             - (exp(-h) + delta1 + delta2) % arma::log1p( exp(h) * AVec )  ) ;
  } else{
    obj_val = arma::accu(  delta1 % ( loglambda01 + eta1 )
                             + (1-delta1) % delta2 % ( loglambda02 + eta2 )
                             + delta1 % delta2 % ( loglambda03 + eta3 ) - AVec );
  }
  return -obj_val;
}



// [[Rcpp::export]]
arma::vec ngradPW_ID(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
                     const arma::vec& delta1, const arma::vec& delta2,
                     const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                     const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
                     const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3,
                     const int frailty_ind){
  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
  int n = X1.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));
  double h;
  if(frailty_ind==1){
    h = para(p01+p02+p03);
  }
  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
  }


  arma::vec Lambda01 = basis1 * arma::exp(phi1);
  arma::vec Lambda02 = basis2 * arma::exp(phi2);
  arma::vec Lambda03 = basis3 * arma::exp(phi3);
  arma::mat BLambda01 = basis1.each_row() % exp(phi1).t();
  arma::mat BLambda02 = basis2.each_row() % exp(phi2).t();
  arma::mat BLambda03 = basis3.each_row() % exp(phi3).t();
  arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);
  arma::vec commonVec;

  arma::vec temp_scorevec(frailty_ind+p01+p02+p03+p1+p2+p3,arma::fill::zeros);
  if(frailty_ind==1){
    commonVec = getCommonVec(delta1, delta2, AVec, h);

    //phi1
    temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * delta1 - BLambda01.t() * (commonVec % arma::exp(h + eta1));
    //phi2
    temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2)
      - BLambda02.t() * (commonVec % arma::exp(h + eta2));
    //phi3
    temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2)
      - BLambda03.t() * (commonVec % arma::exp(h + eta3));
    //h
    temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
    }
    //beta3 (what ina calls u4)
    if(p3 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
        delta1 % commonVec % Lambda03 % arma::exp(h + eta3)); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1)
    }

  } else { //NON-FRAILTY

    //phi1
    temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * delta1 - BLambda01.t() * (arma::exp(eta1));
    //phi2
    temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2)
      - BLambda02.t() * arma::exp(eta2);
    //phi3
    temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2)
      - BLambda03.t() * arma::exp(eta3);
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03, p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - Lambda01 % arma::exp(eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03 + p1, p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - Lambda02 % arma::exp(eta2)) ;
    }
    //beta3 (what ina calls u4)
    if(p3 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03 + p1 + p2,p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
        delta1 % Lambda03 % arma::exp(eta3)); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1)
    }

  }

  return -temp_scorevec;

}




/*

// [[Rcpp::export]]
arma::vec ngradPW_ID_frail(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
                              const arma::vec& delta1, const arma::vec& delta2,
                              const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                              const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
                              const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3){
  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
  int n = X1.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));
  double h = para(p01+p02+p03);

  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
  }

  arma::vec Lambda01 = basis1 * arma::exp(phi1);
  arma::vec Lambda02 = basis2 * arma::exp(phi2);
  arma::vec Lambda03 = basis3 * arma::exp(phi3);
  arma::mat BLambda01 = basis1.each_row() % exp(phi1).t();
  arma::mat BLambda02 = basis2.each_row() % exp(phi2).t();
  arma::mat BLambda03 = basis3.each_row() % exp(phi3).t();
  arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);
  arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

  arma::vec temp_scorevec(1+p01+p02+p03+p1+p2+p3,arma::fill::zeros);
  //phi1
  temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * delta1 - BLambda01.t() * (commonVec % arma::exp(h + eta1));
  //phi2
  temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2)
    - BLambda02.t() * (commonVec % arma::exp(h + eta2));
  //phi3
  temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2)
    - BLambda03.t() * (commonVec % arma::exp(h + eta3));
  //h
  temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));
  //beta1 (what ina calls u2)
  if(p1 > 0){
    temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
  }
  //beta2 (what ina calls u3)
  if(p2 > 0){
    temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
  }
  //beta3 (what ina calls u4)
  if(p3 > 0){
    temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
      delta1 % commonVec % Lambda03 % arma::exp(h + eta3)); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1)
  }
  return -temp_scorevec;

}

*/






/*******
Royston-Parmar
*******/

// [[Rcpp::export]]
double nlogLikRP_uni(const arma::vec& para,
          const arma::vec& y, const arma::vec& delta, const arma::mat& X,
					const arma::mat& basis, const arma::mat& dbasis,
					const arma::mat& basis_yL, const int anyLT){

	//define constants
	int p0 = basis.n_cols;
	int p1 = X.n_cols;
	int n = X.n_rows;

	arma::vec phi = para(arma::span(0, p0-1));

	//define linear predictors
	arma::vec eta(n,arma::fill::zeros);
	if(p1 > 0){
		eta = X * para(arma::span(p0, p0+p1-1));
	}

	//define splines
	arma::vec s = basis * phi;
	arma::vec sprime = dbasis * phi;

	//ll = delta * (log lambda) + log(survival)
	double obj_val; arma::vec exps_yL;
	if(anyLT==1){
	  exps_yL = arma::exp(basis_yL * phi);
	  //zeros in basis_yL * phi indicate need for 0 on exponentiated scale of Lambda0 too
	  exps_yL.replace(1,0);
	  obj_val = arma::accu( delta % (arma::log(sprime) + s + eta - arma::log(y))
                           - (arma::exp(eta) % (arma::exp(s) - exps_yL)) );
	//note that in theory the '-delta*log(y)' term can be proportioned out,
	//but I added it in to keep a properly comparable -2LL scale for AIC
	} else {
    obj_val = arma::accu( delta % (arma::log(sprime) + s + eta - arma::log(y))
                            - arma::exp(s + eta) );
	}

	return(-obj_val);
}

// [[Rcpp::export]]
arma::vec ngradRP_uni(const arma::vec& para,
                      const arma::vec& y, const arma::vec& delta, const arma::mat& X,
                      const arma::mat& basis, const arma::mat& dbasis,
                      const arma::mat& basis_yL, const int anyLT){

	//define constants
	int p0 = basis.n_cols;
	int p1 = X.n_cols;
	int n = X.n_rows;

	arma::vec phi = para(arma::span(0, p0-1));
	//define linear predictors
	arma::vec eta(n, arma::fill::zeros);
	if(p1 > 0){
		eta = X * para(arma::span(p0, p0+p1-1));
	}
	//define splines
	arma::vec s = basis * phi;
	arma::vec sprime = dbasis * phi;

	arma::vec temp_scorevec(p0+p1, arma::fill::zeros);
	arma::vec exps_yL;
	if(anyLT == 1){ //left truncation
	  exps_yL = arma::exp(basis_yL * phi);
	  exps_yL.replace(1,0);
	  //phi
	  temp_scorevec(arma::span(0, p0-1)) = dbasis.t() * (delta / sprime)
                                    	    + basis.t() * (delta - arma::exp(s + eta) )
                                    	    + basis_yL.t() * (exps_yL % arma::exp(eta));
    //beta
    if(p1 > 0){
      temp_scorevec(arma::span(p0,p0+p1-1)) = X.t() * ( delta - arma::exp(eta) % (arma::exp(s) - exps_yL) );
    }
	} else { //no left truncation
	  //phi
	  temp_scorevec(arma::span(0, p0-1)) = dbasis.t() * (delta / sprime)
                                    	    + basis.t() * (delta - arma::exp(s + eta) );
	  //beta
	  if(p1 > 0){
	    temp_scorevec(arma::span(p0,p0+p1-1)) = X.t() * ( delta - arma::exp(s + eta) );
	  }
	}
	return(-temp_scorevec);
}

// [[Rcpp::export]]
double nlogLikRP_ID(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
                             const arma::vec& delta1, const arma::vec& delta2,
                             const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                             const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, const arma::mat& basis3_y1,
                             const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3,
                             const std::string model, const int frailty_ind){

  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
  int n = X1.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));

  double h;
  if(frailty_ind==1){
    h = para(p01+p02+p03);
  }
  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
  }



  //define splines
  //assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
  arma::vec s1 = basis1 * phi1;
  arma::vec s1prime = dbasis1 * phi1;
  //assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
  arma::vec s2 = basis2 * phi2;
  arma::vec s2prime = dbasis2 * phi2;
  //assumptions: under semi-markov basis 3 is coming from y_2-y1, with placement of knots depending on quantiles of (y2-y1)[delta1==1 & delta2==1]
  arma::vec s3 = basis3 * phi3;
  arma::vec s3prime = delta1 % (dbasis3 * phi3);
  s3prime.replace(0,1);

  arma::vec AVec,s3_y1,logdiff;
  if (model.compare("markov") == 0){ //markov
    //if we want a markov approach, we also need extra piece, which is y_1 set at same knots of y2[delta1==1 & delta2==1]
    arma::vec s3_y1 = basis3_y1 * phi3;
    AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % (arma::exp(s3 + eta3) - arma::exp(s3_y1 + eta3));
    logdiff = arma::log(y2); //awkward name reflects the fact that the markov assumption sets log(y2|y1) = log(y2), so there's no 'difference' here
  } else{
    AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % arma::exp(s3 + eta3);
    //you might think that all of the information about y1 and y2 is coming through the basis, but
    //to keep the likelihood on the same scale as other parametric models, we need to include the 'extraneous'
    //terms involving y1, y2, and y2-y1 in the final calculation
    logdiff = arma::log(y2-y1);
    logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.
  }

  double obj_val;
  if(frailty_ind==1){
    obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1))
                             + (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
                             + delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - logdiff + log1p( exp(h) )  )
                             - (exp(-h) + delta1 + delta2) % arma::log1p( exp(h) * AVec )  ) ;
  } else {
    obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1))
                             + (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
                             + delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - logdiff  ) - AVec );
  }
  return -obj_val;
}


/*

// [[Rcpp::export]]
double nlogLikRP_ID_frail_SM(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
							 const arma::vec& delta1, const arma::vec& delta2,
						   	 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
						   	 const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
						   	 const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3,
						   	 const int frailty_ind){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1));
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));

	double h;
	if(frailty_ind==1){
	  h = para(p01+p02+p03);
	}
	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
	  eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
	}
	if(p2 > 0){
	  eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
	}
	if(p3 > 0){
	  eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
	}

	//you might think that all of the information about y1 and y2 is coming through the basis, but
	//to keep the likelihood on the same scale as other parametric models, we need to include the 'extraneous'
	//terms involving y1, y2, and y2-y1 in the final calculation
	arma::vec logdiff = arma::log(y2-y1);
	logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.


	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under semi-markov basis 3 is coming from y_2-y1, with placement of knots depending on quantiles of (y2-y1)[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = delta1 % (dbasis3 * phi3);
	s3prime.replace(0,1);

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2,
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec;
	AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % arma::exp(s3 + eta3);

	double obj_val;
	if(frailty_ind==1){
	  obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1))
                            + (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
                            + delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - logdiff + log1p( exp(h) )  )
                            - (exp(-h) + delta1 + delta2) % arma::log1p( exp(h) * AVec )  ) ;
	} else {
	  obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1))
                            + (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
                            + delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - logdiff  ) - AVec );
	}

  return -obj_val;
}


// [[Rcpp::export]]
double nlogLikRP_ID_frail_M(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
							const arma::vec& delta1, const arma::vec& delta2,
						    const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
						    const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, const arma::mat& basis3_y1,
						    const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3,
						    const int frailty_ind){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1));
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));

	double h;
	if(frailty_ind==1){
	  h = para(p01+p02+p03);
	}
	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
	  eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
	}
	if(p2 > 0){
	  eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
	}
	if(p3 > 0){
	  eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
	}

	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under markov, basis 3 is coming from y_2, with placement of knots depending on quantiles of y2[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = dbasis3 * phi3;

	//if we want a markov approach, we also need extra piece, which is y_1 set at same knots of y2[delta1==1 & delta2==1]
	arma::vec s3_y1 = basis3_y1 * phi3;

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2,
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec;
	AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % (arma::exp(s3 + eta3) - arma::exp(s3_y1 + eta3));

	double obj_val;
	if(frailty_ind == 1){
	  obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1) )
                            + (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
                            + delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - arma::log(y2) + log1p( exp(h) )  )
                            - (exp(-h) + delta1 + delta2) % arma::log1p( exp(h) * AVec )  ) ;
	} else {
	  obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1) )
                            + (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
                            + delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - arma::log(y2)  ) - AVec );
	}

  return -obj_val;
}

*/



// [[Rcpp::export]]
arma::vec ngradRP_ID(const arma::vec& para, const arma::vec& y1, const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
                    const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                    const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, const arma::mat& basis3_y1,
                    const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3,
                    const std::string model, const int frailty_ind){
  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
  int n = X1.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));
  double h;
  if(frailty_ind==1){
    h = para(p01+p02+p03);
  }
  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
  }

  //define splines
  //assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
  arma::vec s1 = basis1 * phi1;
  arma::vec s1prime = dbasis1 * phi1;
  //assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
  arma::vec s2 = basis2 * phi2;
  arma::vec s2prime = dbasis2 * phi2;
  //assumptions: under semi-markov basis 3 is coming from y_2-y1, with placement of knots depending on quantiles of (y2-y1)[delta1==1 & delta2==1]
  arma::vec s3 = basis3 * phi3;
  arma::vec s3prime = delta1 % (dbasis3 * phi3);
  s3prime.replace(0,1);

  arma::vec AVec,s3_y1;
  if (model.compare("markov") == 0){ //markov
    //if we want a markov approach, we also need extra piece, which is y_1 set at same knots of y2[delta1==1 & delta2==1]
    s3_y1 = basis3_y1 * phi3;
    AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % (arma::exp(s3 + eta3) - arma::exp(s3_y1 + eta3));
  } else{
    AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % arma::exp(s3 + eta3);
  }

  arma::vec commonVec;
  arma::vec temp_scorevec(frailty_ind+p01+p02+p03+p1+p2+p3,arma::fill::zeros);
  if(frailty_ind==1){
    commonVec = getCommonVec(delta1, delta2, AVec, h);

    //phi1
    temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * (delta1 / s1prime)
      + basis1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1) );
    //phi2
    temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2 / s2prime)
      + basis2.t() * ((1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2) );
    //h
    temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2)) ;
    }

    if (model.compare("markov") == 0){ //markov
      //phi3
      temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 / s3prime)
      + basis3.t() * (delta1 % delta2)
      - (basis3.t() * (delta1 % commonVec % arma::exp(h + s3 + eta3))
      - basis3_y1.t() * (delta1 % commonVec % arma::exp(h + s3_y1 + eta3)));  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.
      //beta3
      if(p3 > 0){
        temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
          delta1 % commonVec % arma::exp(h + eta3) % (arma::exp(s3) - arma::exp(s3_y1))); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.
      }
    } else {
      //phi3
      temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 / s3prime)
      + basis3.t() * (delta1 % delta2 -
        delta1 % commonVec % arma::exp(h + s3 + eta3) );  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1). Helps make it more robust
      //beta3 (what ina calls u4)
      if(p3 > 0){
        temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
          delta1 % commonVec % arma::exp(s3 + h + eta3)); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1)
      }
    }

  } else { //NON-FRAILTY
    //phi1
    temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * (delta1 / s1prime)
      + basis1.t() * (delta1 - arma::exp(s1 + eta1) );
    //phi2
    temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2 / s2prime)
      + basis2.t() * ((1-delta1) % delta2 - arma::exp(s2 + eta2) );
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03, p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - arma::exp(s1 + eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03 + p1, p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - arma::exp(s2 + eta2)) ;
    }

    if (model.compare("markov") == 0){ //markov
      //phi3
      temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 / s3prime)
      + basis3.t() * (delta1 % delta2)
      - (basis3.t() * (delta1 % arma::exp(s3 + eta3))
      - basis3_y1.t() * (delta1 % arma::exp(s3_y1 + eta3)));  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.
      //beta3
      if(p3 > 0){
        temp_scorevec(arma::span(p01 + p02 + p03 + p1 + p2, p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
          delta1 % arma::exp(eta3) % (arma::exp(s3) - arma::exp(s3_y1))); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.
      }

    } else { //semi-markov
      //phi3
      temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 / s3prime)
      + basis3.t() * (delta1 % delta2 -
        delta1 % arma::exp(s3 + eta3) );  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1). Helps make it more robust
      //beta3 (what ina calls u4)
      if(p3 > 0){
        temp_scorevec(arma::span(p01 + p02 + p03 + p1 + p2,p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
          delta1 % arma::exp(s3 + eta3)); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1)
      }
    }

  }
  return -temp_scorevec;
}


/*


// [[Rcpp::export]]
arma::vec ngradRP_ID_frail_SM(const arma::vec& para, const arma::vec& y1, const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
						   	  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
						      const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
						      const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1));
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));
	double h = para(p01+p02+p03);

	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
		eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
	}
	if(p2 > 0){
		eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
	}
	if(p3 > 0){
		eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
	}

	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under semi-markov basis 3 is coming from y_2-y1, with placement of knots depending on quantiles of (y2-y1)[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = dbasis3 * phi3;

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2,
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % arma::exp(s3 + eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

	arma::vec temp_scorevec(1+p01+p02+p03+p1+p2+p3,arma::fill::zeros);
	//phi1
	temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * (delta1 / s1prime)
											+ basis1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1) );
	//phi2
	temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2 / s2prime)
													+ basis2.t() * ((1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2) );
	//phi3
	temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 / s3prime)
													+ basis3.t() * (delta1 % delta2 -
														delta1 % commonVec % arma::exp(h + s3 + eta3) );  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1). Helps make it more robust
	//h
	temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));
	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1));
	}
	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2)) ;
	}
	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
								 delta1 % commonVec % arma::exp(s3 + h + eta3)); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1)
	}
  return -temp_scorevec;
}



//This function implements both markov and semi-markov specifications, but I think I'll need to write separate functions for the grad and hess
// [[Rcpp::export]]
arma::vec ngradRP_ID_frail_M(const arma::vec& para,
                 const arma::vec& y1, const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
						     const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
						     const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, const arma::mat& basis3_y1,
						     const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1));
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));
	double h = para(p01+p02+p03);

	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
		eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
	}
	if(p2 > 0){
		eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
	}
	if(p3 > 0){
		eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
	}

	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under markov, basis 3 is coming from y_2, with placement of knots depending on quantiles of y2[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = dbasis3 * phi3;
	//if we want a markov approach, we also need extra piece, which is y_1 set at same knots of y2[delta1==1 & delta2==1]
	arma::vec s3_y1 = basis3_y1 * phi3;

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2,
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % (arma::exp(s3 + eta3) - arma::exp(s3_y1 + eta3));
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

  arma::vec temp_scorevec(1+p01+p02+p03+p1+p2+p3,arma::fill::zeros);
	//phi1
	temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * (delta1 / s1prime)
											+ basis1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1) );
	//phi2
	temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2 / s2prime)
													+ basis2.t() * ((1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2) );
	//phi3
	temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 / s3prime)
													+ basis3.t() * (delta1 % delta2)
													- (basis3.t() * (delta1 % commonVec % arma::exp(h + s3 + eta3))
													   - basis3_y1.t() * (delta1 % commonVec % arma::exp(h + s3_y1 + eta3)));  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.
	//h
	temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));
	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1));
	}
	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2)) ;
	}
	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 -
								 delta1 % commonVec % arma::exp(h + eta3) % (arma::exp(s3) - arma::exp(s3_y1))); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.
	}
  return -temp_scorevec;
}

 */

/*******
 B-Spline
 *******/

arma::vec getLambda0BS(const arma::vec& y,
                       const arma::vec& lambda0,const arma::vec& quad_weights){
  int n = y.n_rows;
  int n_quad = quad_weights.n_rows;
  arma::vec temp(n,arma::fill::zeros);
  for(int i = 0; i < n; i++) {
    if(y(i)>0){
      temp(i) = arma::dot( lambda0(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)),quad_weights) * y(i) / 2;
    }
  }
  return temp;
}

arma::vec getAVecBS(const arma::vec& y1,const arma::vec& y2,
                    const arma::vec& delta1,const arma::vec& delta2,
                    const arma::vec& lambda01, const arma::vec& lambda02, const arma::vec& lambda03,
                    const arma::vec& eta1, const arma::vec& eta2, const arma::vec& eta3,
                    const arma::vec& quad_weights){

  int n = y1.n_rows;
  int n_quad = quad_weights.n_rows;
  arma::vec temp(n);

  //loglambda vectors are of length (n_quad + 1) * n, with the first n corresponding to
  //the actual event time y1/y2/y2-y1, and then each subsequent set of 15 corresponding
  //to
  for(int i = 0; i < n; i++) {
    temp(i) = arma::dot( lambda01(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)),quad_weights)
    * exp(eta1(i)) * y1(i) / 2
    + arma::dot( lambda02(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)),quad_weights)
    * exp(eta2(i)) * y1(i) / 2;
    if(delta1(i)==1){
      temp(i) += arma::dot( lambda03(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)),quad_weights)
      * exp(eta3(i)) * (y2(i) - y1(i)) / 2;
    }
  }
  return temp;
}

// [[Rcpp::export]]
double nlogLikBS_uni(const arma::vec& para,
                     const arma::vec& y, const arma::vec& delta, const arma::mat& X,
                     const arma::mat& basis, const arma::vec& quad_weights){

  //define constants
  int p0 = basis.n_cols;
  int p1 = X.n_cols;
  int n = X.n_rows;

  arma::vec phi = para(arma::span(0, p0-1));
  arma::vec eta(n,arma::fill::zeros);
  if(p1 > 0){
    eta = X * para(arma::span(p0, p0+p1-1));
  }
  arma::vec loglambda0 = basis * phi;

  //delta * (log lambda) + log(survival)
  double obj_val = arma::accu( delta % (loglambda0(arma::span(0,n-1)) + eta)
                      - getLambda0BS(y,arma::exp(loglambda0),quad_weights)
                        % arma::exp(eta) );

  return(-obj_val);
}

// [[Rcpp::export]]
arma::vec ngradBS_uni(const arma::vec& para,
                   const arma::vec& y, const arma::vec& delta, const arma::mat& X,
                   const arma::mat& basis, const arma::vec& quad_weights){

  //define constants
  int p0 = basis.n_cols;
  int p1 = X.n_cols;
  int n = X.n_rows;
  int n_quad = quad_weights.n_rows;

  arma::vec phi = para(arma::span(0, p0-1));
  arma::vec eta(n,arma::fill::zeros);
  if(p1 > 0){
    eta = X * para(arma::span(p0, p0+p1-1));
  }
  arma::vec lambda0 = arma::exp(basis * phi);
  arma::vec Lambda0 = getLambda0BS(y,lambda0,quad_weights);

  //matrix with each column representing the integralint_0^{t_i} B(s) * h_0(s) ds, and each row is an i
  arma::mat BLambda0(n,p0,arma::fill::zeros);
  arma::mat basisi(n_quad,p0);
  for(int i = 0; i < n; i++) {
    basisi = basis.rows(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1));
    //multiply each basis column by lambda, yielding an n_quad by k matrix. Then multiply with 1 by n_quad weight vector
    BLambda0.row(i) = y(i)/2 * quad_weights.t() *
      ( basisi.each_col()
          % lambda0(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)) );
  }

  arma::vec temp_scorevec(p0+p1, arma::fill::zeros);
  temp_scorevec(arma::span(0, p0-1)) = basis.rows(0,n-1).t() * delta - BLambda0.t() * arma::exp(eta);
  if(p1 > 0){
    temp_scorevec(arma::span(p0,p0+p1-1)) = X.t() * (delta - Lambda0 % arma::exp(eta));
  }
  return -temp_scorevec;
}

// [[Rcpp::export]]
double nlogLikBS_ID(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
                    const arma::vec& delta1, const arma::vec& delta2,
                    const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                    const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
                    const arma::vec& quad_weights, const int frailty_ind){
  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols; int n = X1.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));

  double h;
  if(frailty_ind==1){
    h = para(p01+p02+p03);
  }

  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
  }

  //first n rows are the observed event times,
  //then every subsequent set of 15 rows corresponds with one subject
  arma::vec loglambda01 = basis1 * phi1;
  arma::vec loglambda02 = basis2 * phi2;
  arma::vec loglambda03 = basis3 * phi3;
  arma::vec AVec = getAVecBS(y1,y2,delta1,delta2,
                             arma::exp(loglambda01),arma::exp(loglambda02),arma::exp(loglambda03),
                             eta1,eta2,eta3,quad_weights);

  double obj_val;
  if(frailty_ind==1){
   obj_val = arma::accu(  delta1 % (loglambda01(arma::span(0,n-1)) + eta1)
                            + (1-delta1) % delta2 % (loglambda02(arma::span(0,n-1)) + eta2)
                            + delta1 % delta2 % (loglambda03(arma::span(0,n-1)) + eta3 + log1p(exp(h)))
                            - (exp(-h) + delta1 + delta2) % arma::log1p( exp(h) * AVec ) );
  } else {
    obj_val = arma::accu(  delta1 % (loglambda01(arma::span(0,n-1)) + eta1)
                             + (1-delta1) % delta2 % (loglambda02(arma::span(0,n-1)) + eta2)
                             + delta1 % delta2 % (loglambda03(arma::span(0,n-1)) + eta3) - AVec );
  }
  return -obj_val;
}



// [[Rcpp::export]]
arma::vec ngradBS_ID(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
                           const arma::vec& delta1, const arma::vec& delta2,
                           const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                           const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
                           const arma::vec& quad_weights, const int frailty_ind){
  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
  int n = X1.n_rows; int n_quad = quad_weights.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));
  double h;
  if(frailty_ind==1){
    h = para(p01+p02+p03);
  }

  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(frailty_ind+p01+p02+p03,frailty_ind-1+p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(frailty_ind+p01+p02+p03+p1, frailty_ind-1+p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(frailty_ind+p01+p02+p03+p1+p2, frailty_ind-1+p01+p02+p03+p1+p2+p3));
  }

  arma::vec lambda01 = exp(basis1 * phi1);
  arma::vec lambda02 = exp(basis2 * phi2);
  arma::vec lambda03 = exp(basis3 * phi3);

  //these vectors are of length n * (n_quad + 1)
  //first n rows are the observed event times,
  //then every subsequent set of 15 rows corresponds with
  //quadrature points for one subject
  arma::vec Lambda01 = getLambda0BS(y1,lambda01,quad_weights);
  arma::vec Lambda02 = getLambda0BS(y1,lambda02,quad_weights);
  arma::vec Lambda03 = getLambda0BS(y2-y1,lambda03,quad_weights);
  //now, generate a length n vector of the sum of cumulative hazards for each subject
  arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);

  //matrix with each column representing the integralint_0^{t_i} B(s) * h_0(s) ds, and each row is an i
  arma::mat BLambda01(n,p01,arma::fill::zeros);
  arma::mat BLambda02(n,p02,arma::fill::zeros);
  arma::mat BLambda03(n,p03,arma::fill::zeros);

  arma::mat basis1i(n_quad,p01);
  arma::mat basis2i(n_quad,p02);
  arma::mat basis3i(n_quad,p03);
  for(int i = 0; i < n; i++) {

    basis1i = basis1.rows(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1));
    basis2i = basis2.rows(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1));
    basis3i = basis3.rows(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1));

    //multiply each basis column by lambda, yielding an n_quad by k matrix. Then multiply with 1 by n_quad weight vector
    BLambda01.row(i) = y1(i)/2 * quad_weights.t() *
      ( basis1i.each_col()
          % lambda01(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)) );
    BLambda02.row(i) = y1(i)/2 * quad_weights.t() *
    ( basis2i.each_col()
        % lambda02(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)) );
    if(delta1(i)==1){
      BLambda03.row(i) = (y2(i)-y1(i))/2 * quad_weights.t() *
        ( basis3i.each_col()
            % lambda03(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)) );
    }
  }

  arma::vec commonVec;
  arma::vec temp_scorevec(frailty_ind+p01+p02+p03+p1+p2+p3,arma::fill::zeros);
  if(frailty_ind==1){
    commonVec = getCommonVec(delta1, delta2, AVec, h);
    //phi1
    temp_scorevec(arma::span(0, p01 - 1)) = basis1.rows(0,n-1).t() * delta1 - BLambda01.t() * (commonVec % arma::exp(h + eta1)) ;
    //phi2
    temp_scorevec(arma::span(p01, p01 + p02 - 1)) = basis2.rows(0,n-1).t() * ((1-delta1) % delta2) - BLambda02.t() * (commonVec % arma::exp(h + eta2));
    //phi3
    temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = basis3.rows(0,n-1).t() * (delta1 % delta2) - BLambda03.t() * (commonVec % arma::exp(h + eta3));
    //h
    temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
    }
    //beta3 (what ina calls u4)
    if(p3 > 0){
      temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
    }
  } else {
    //phi1
    temp_scorevec(arma::span(0, p01 - 1)) = basis1.rows(0,n-1).t() * delta1 - BLambda01.t() * arma::exp(eta1) ;
    //phi2
    temp_scorevec(arma::span(p01, p01 + p02 - 1)) = basis2.rows(0,n-1).t() * ((1-delta1) % delta2) - BLambda02.t() * arma::exp(eta2);
    //phi3
    temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = basis3.rows(0,n-1).t() * (delta1 % delta2) - BLambda03.t() * arma::exp(eta3);
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03, p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - Lambda01 % arma::exp(eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03 + p1, p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - Lambda02 % arma::exp(eta2)) ;
    }
    //beta3 (what ina calls u4)
    if(p3 > 0){
      temp_scorevec(arma::span(p01 + p02 + p03 + p1 + p2,p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - Lambda03 % arma::exp(eta3));
    }
  }

  return -temp_scorevec;
}


/*

// [[Rcpp::export]]
arma::vec ngradBS_ID_frail(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
                           const arma::vec& delta1, const arma::vec& delta2,
                           const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                           const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3,
                           const arma::vec& quad_weights){
  //define constants
  int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;
  int n = X1.n_rows; int n_quad = quad_weights.n_rows;

  arma::vec phi1 = para(arma::span(0, p01-1));
  arma::vec phi2 = para(arma::span(p01, p01+p02-1));
  arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1));
  double h = para(p01+p02+p03);

  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
  }

  arma::vec lambda01 = exp(basis1 * phi1);
  arma::vec lambda02 = exp(basis2 * phi2);
  arma::vec lambda03 = exp(basis3 * phi3);

  //these vectors are of length n * (n_quad + 1)
  //first n rows are the observed event times,
  //then every subsequent set of 15 rows corresponds with
  //quadrature points for one subject
  arma::vec Lambda01 = getLambda0BS(y1,lambda01,quad_weights);
  arma::vec Lambda02 = getLambda0BS(y1,lambda02,quad_weights);
  arma::vec Lambda03 = getLambda0BS(y2-y1,lambda03,quad_weights);
  //now, generate a length n vector of the sum of cumulative hazards for each subject
  arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);
  arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

  //matrix with each column representing the integralint_0^{t_i} B(s) * h_0(s) ds, and each row is an i
  arma::mat BLambda01(n,p01,arma::fill::zeros);
  arma::mat BLambda02(n,p02,arma::fill::zeros);
  arma::mat BLambda03(n,p03,arma::fill::zeros);

  arma::mat basis1i(n_quad,p01);
  arma::mat basis2i(n_quad,p02);
  arma::mat basis3i(n_quad,p03);
  for(int i = 0; i < n; i++) {

    basis1i = basis1.rows(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1));
    basis2i = basis2.rows(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1));
    basis3i = basis3.rows(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1));

    //multiply each basis column by lambda, yielding an n_quad by k matrix. Then multiply with 1 by n_quad weight vector
    BLambda01.row(i) = y1(i)/2 * quad_weights.t() *
      ( basis1i.each_col()
          % lambda01(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)) );
    BLambda02.row(i) = y1(i)/2 * quad_weights.t() *
    ( basis2i.each_col()
        % lambda02(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)) );
    if(delta1(i)==1){
      BLambda03.row(i) = (y2(i)-y1(i))/2 * quad_weights.t() *
        ( basis3i.each_col()
            % lambda03(arma::span(n + i*n_quad,n + (i+1)*n_quad - 1)) );
    }
  }

  arma::vec temp_scorevec(1+p01+p02+p03+p1+p2+p3,arma::fill::zeros);
  //phi1
  temp_scorevec(arma::span(0, p01 - 1)) = basis1.rows(0,n-1).t() * delta1 - BLambda01.t() * (commonVec % arma::exp(h + eta1)) ;
  //phi2
  temp_scorevec(arma::span(p01, p01 + p02 - 1)) = basis2.rows(0,n-1).t() * ((1-delta1) % delta2) - BLambda02.t() * (commonVec % arma::exp(h + eta2));
  //phi3
  temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = basis3.rows(0,n-1).t() * (delta1 % delta2) - BLambda03.t() * (commonVec % arma::exp(h + eta3));
  //h
  temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));
  //beta1 (what ina calls u2)
  if(p1 > 0){
    temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
  }
  //beta2 (what ina calls u3)
  if(p2 > 0){
    temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
  }
  //beta3 (what ina calls u4)
  if(p3 > 0){
    temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
  }
  return -temp_scorevec;
}


*/



/*******
WEIBULL
*******/


arma::vec getLambda0WB(const arma::vec& y, double a, double k){
	return exp(k) * arma::pow(y, exp(a));
}

// [[Rcpp::export]]
double nlogLikWB_uni(const arma::vec& para, const arma::vec& y,
                     const arma::vec& delta, const arma::mat& X,
                     const arma::vec& yL, const int anyLT){

  //define constants
  int p1 = X.n_cols; int n = X.n_rows;
  double k = para(0); double a = para(1);
  arma::vec eta(n,arma::fill::zeros);
  if(p1 > 0){
    eta = X * para(arma::span(2, 2+p1-1));
  }
  //delta * (log lambda) + log(survival)
  double obj_val = arma::accu( delta % (a + k + (exp(a) - 1) * arma::log(y) + eta)
                                 - getLambda0WB(y,a,k) % arma::exp(eta) );
  if(anyLT==1){
    obj_val += arma::accu( getLambda0WB(yL,a,k) % arma::exp(eta) );
  }
  return(-obj_val);
}

// [[Rcpp::export]]
arma::vec ngradWB_uni(const arma::vec& para, const arma::vec& y,
                     const arma::vec& delta, const arma::mat& X,
                     const arma::vec& yL, const int anyLT){

  //define constants
  int p1 = X.n_cols; int n = X.n_rows;
  double k = para(0); double a = para(1);
  arma::vec eta(n,arma::fill::zeros);
  if(p1 > 0){
    eta = X * para(arma::span(2, 2+p1-1));
  }
  arma::vec Lambda0 = getLambda0WB(y,a,k);

  arma::vec temp_scorevec(2+p1, arma::fill::zeros);
  arma::vec Lambda0L, logyL; //quantities for left truncation
  if(anyLT==1){ //left truncation
    Lambda0L = getLambda0WB(yL,a,k);
    logyL = arma::log(yL);
    logyL = logyL.replace(-arma::datum::inf, 0); //log(yL) with the negative infinity values replaced with 0's.

    temp_scorevec(0) = arma::accu( delta - (Lambda0-Lambda0L) % arma::exp(eta) );
    temp_scorevec(1) = arma::accu( delta % (1 + exp(a) * arma::log(y)) -
      arma::exp(a+eta) % (Lambda0 % arma::log(y) - Lambda0L % logyL) );
    if(p1 > 0){
      temp_scorevec(arma::span(2,2+p1-1)) = X.t() * (delta - (Lambda0-Lambda0L) % arma::exp(eta));
    }
  } else{
    temp_scorevec(0) = arma::accu( delta - Lambda0 % arma::exp(eta) );
    temp_scorevec(1) = arma::accu( delta % (1 + exp(a) * arma::log(y)) -
      Lambda0 % arma::exp(a+eta) % arma::log(y));
    if(p1 > 0){
      temp_scorevec(arma::span(2,2+p1-1)) = X.t() * (delta - Lambda0 % arma::exp(eta));
    }
  }

  return(-temp_scorevec);
}






// [[Rcpp::export]]
double nlogLikWB_ID(const arma::vec& para,
						 const arma::vec& y1,const arma::vec& y2,
						 const arma::vec& delta1, const arma::vec& delta2,
						 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
						 const std::string model, const int frailty_ind){
  //This function implements both markov and semi-markov specifications,
  //but I think I'll need to write separate functions for the grad and hess

  //define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);

	double h;
	if(frailty_ind==1){
	  h = para(6);
	}

	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
		eta1 = X1 * para(arma::span(frailty_ind + 6, frailty_ind + 5 + p1));
	}
	if(p2 > 0){
		eta2 = X2 * para(arma::span(frailty_ind + 6 + p1, frailty_ind + 5 + p1 + p2));
	}
	if(p3 > 0){
		eta3 = X3 * para(arma::span(frailty_ind + 6 + p1 + p2 , frailty_ind + 5 + p1 + p2 + p3));
	}

	//define collected terms
	arma::vec AVec, logdiff;
	if (model.compare("markov") == 0){ //markov
		AVec = getLambda0WB(y1,a1,k1) % arma::exp(eta1)
	       + getLambda0WB(y1,a2,k2) % arma::exp(eta2)
	       + (getLambda0WB(y2,a3,k3) - getLambda0WB(y1,a3,k3)) % arma::exp(eta3);
		logdiff = arma::log(y2); //awkward name reflects the fact that the markov assumption sets log(y2|y1) = log(y2), so there's no 'difference' here
	} else { //semi-markov
		AVec = getLambda0WB(y1,a1,k1) % arma::exp(eta1)
	       + getLambda0WB(y1,a2,k2) % arma::exp(eta2)
	       + getLambda0WB(y2-y1,a3,k3) % arma::exp(eta3);
		logdiff = arma::log(y2-y1);
		logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.
	}


	double obj_val;
	if(frailty_ind ==1){
	  obj_val = arma::accu( delta1 % ( a1 + k1 + (exp(a1) - 1) * arma::log(y1) + eta1 )
                           + (1-delta1) % delta2 %  (a2 + k2 + (exp(a2) - 1) * arma::log(y1) + eta2)
                           + delta1 % delta2 % (a3 + k3 + (exp(a3) - 1) * logdiff + eta3 + log1p(exp(h)))
                           - (exp(-h) + delta1 + delta2) % arma::log1p(exp(h) * AVec) );
	} else {
	  obj_val = arma::accu( delta1 % ( a1 + k1 + (exp(a1) - 1) * arma::log(y1) + eta1 )
                           + (1-delta1) % delta2 %  (a2 + k2 + (exp(a2) - 1) * arma::log(y1) + eta2)
                           + delta1 % delta2 % (a3 + k3 + (exp(a3) - 1) * logdiff + eta3) - AVec );
	}

	// //this is the slow way of writing it, I could simplify to make this a bit quicker
	// arma::vec LL1 = delta1 % delta2 % (  a1 + k1 + (exp(a1) - 1) * arma::log(y1) + eta1
	// 								   + a3 + k3 + (exp(a3) - 1) * logdiff       + eta3
	// 								   + log1p(exp(h)) //added because I think it's missing
	// 								   - (exp(-h) + 2) * arma::log1p(exp(h) * AVec)
	// 								  ); //LL1
	// arma::vec LL3 = delta1 % (1-delta2) % (a1 + k1 + (exp(a1) - 1) * arma::log(y1) + eta1
	// 										- (exp(-h)+1) * arma::log1p(exp(h) * AVec)); //LL3	;
	// arma::vec LL2 = (1-delta1) % delta2 % (a2 + k2 + (exp(a2) - 1) * arma::log(y1) + eta2
	// 										- (exp(-h)+1) * arma::log1p(exp(h) * AVec)); // LL 2
	// arma::vec LL4 = (1-delta1) % (1-delta2) % (-exp(-h) * arma::log1p(exp(h) * AVec)); // LL4
	// double obj_val = arma::accu( LL1 + LL2 + LL3 + LL4);

  return -obj_val;
}


// [[Rcpp::export]]
arma::vec ngradWB_ID(const arma::vec& para,
                    const arma::vec& y1,const arma::vec& y2,
                    const arma::vec& delta1, const arma::vec& delta2,
                    const arma::mat& X1, const arma::mat& X2, const arma::mat& X3,
                    const std::string model, const int frailty_ind){
  //define constants
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

  double k1 = para(0);
  double a1 = para(1);
  double k2 = para(2);
  double a2 = para(3);
  double k3 = para(4);
  double a3 = para(5);
  double h;
  if(frailty_ind==1){
    h = para(6);
  }
  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(frailty_ind + 6, frailty_ind + 5 + p1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(frailty_ind + 6 + p1, frailty_ind + 5 + p1 + p2));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(frailty_ind + 6 + p1 + p2 , frailty_ind + 5 + p1 + p2 + p3));
  }

  arma::vec Lambda01 = getLambda0WB(y1, a1, k1);
  arma::vec Lambda02 = getLambda0WB(y1, a2, k2);
  arma::vec AVec, logdiff, Lambda03,Lambda03y2,Lambda03y1, commonVec;
  if (model.compare("markov") == 0){ //markov
    Lambda03y2 = getLambda0WB(y2, a3, k3);
    Lambda03y1 = getLambda0WB(y1, a3, k3);
    AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2)
      + (Lambda03y2-Lambda03y1) % arma::exp(eta3);
  } else { //semi-markov
    Lambda03 = getLambda0WB(y2-y1, a3, k3);
    AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);
    logdiff = arma::log(y2-y1);
    logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.
  }

  arma::vec temp_scorevec(p1+p2+p3+6+frailty_ind,arma::fill::zeros);
  if(frailty_ind==1){
    commonVec = getCommonVec(delta1, delta2, AVec, h);
    //k1 (what ina calls u5)
    temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
    //a1 (what ina calls u6)
    temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) -
      commonVec % Lambda01 % arma::exp(h+a1+eta1) % arma::log(y1));
    //k2 (what ina calls u7)
    temp_scorevec(2) = arma::accu( (1-delta1) % delta2 -
      commonVec % Lambda02 % arma::exp(h+eta2));
    //a2 (what ina calls u8)
    temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) -
      commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );
    //h (what ina calls u1)
    temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec));
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
    }

    if (model.compare("markov") == 0){ //markov
      //k3 (what ina calls u9)
      temp_scorevec(4) = arma::accu( delta1 % delta2 -
        commonVec % (Lambda03y2 - Lambda03y1) % arma::exp(h+eta3));
      //a3 (what ina calls u10)
      temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) -
        commonVec % arma::exp(h+a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) );
      //beta3 (what ina calls u4)
      if(p3 > 0){
        temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % (Lambda03y2 - Lambda03y1)  % arma::exp(h + eta3));
      }
    } else {
      //k3 (what ina calls u9)
      temp_scorevec(4) = arma::accu( delta1 % delta2 -
        commonVec % Lambda03 % arma::exp(h+eta3));
      //a3 (what ina calls u10)
      temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * logdiff) -
        commonVec % Lambda03 % arma::exp(h+a3+eta3) % logdiff );
      //beta3 (what ina calls u4)
      if(p3 > 0){
        temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
      }
    }

  } else { //NON-FRAILTY
    //k1 (what ina calls u5)
    temp_scorevec(0) = arma::accu(delta1 - Lambda01 % arma::exp(eta1));
    //a1 (what ina calls u6)
    temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) -
      Lambda01 % arma::exp(a1+eta1) % arma::log(y1));
    //k2 (what ina calls u7)
    temp_scorevec(2) = arma::accu( (1-delta1) % delta2 - Lambda02 % arma::exp(eta2));
    //a2 (what ina calls u8)
    temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) -
      Lambda02 % arma::exp(a2+eta2) % arma::log(y1) );
    //beta1 (what ina calls u2)
    if(p1 > 0){
      temp_scorevec(arma::span(6, 6 + p1 - 1)) = X1.t() * (delta1 - Lambda01 % arma::exp(eta1));
    }
    //beta2 (what ina calls u3)
    if(p2 > 0){
      temp_scorevec(arma::span(6 + p1, 6 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - Lambda02 % arma::exp(eta2)) ;
    }

    if (model.compare("markov") == 0){ //markov
      //k3 (what ina calls u9)
      temp_scorevec(4) = arma::accu( delta1 % delta2 -
        (Lambda03y2 - Lambda03y1) % arma::exp(eta3));
      //a3 (what ina calls u10)
      temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) -
        arma::exp(a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) );
      //beta3 (what ina calls u4)
      if(p3 > 0){
        temp_scorevec(arma::span(6 + p1 + p2 , 6 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - (Lambda03y2 - Lambda03y1)  % arma::exp(eta3));
      }
    } else {
      //k3 (what ina calls u9)
      temp_scorevec(4) = arma::accu( delta1 % delta2 -
        Lambda03 % arma::exp(eta3));
      //a3 (what ina calls u10)
      temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * logdiff) -
        Lambda03 % arma::exp(a3+eta3) % logdiff );
      //beta3 (what ina calls u4)
      if(p3 > 0){
        temp_scorevec(arma::span(6 + p1 + p2 , 6 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - Lambda03 % arma::exp(eta3));
      }
    }
  }

  return -temp_scorevec;
}


/*

// [[Rcpp::export]]
arma::vec ngradWB_ID_frail_SM(const arma::vec& para,
							  const arma::vec& y1,const arma::vec& y2,
							  const arma::vec& delta1, const arma::vec& delta2,
							  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
	  eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 > 0){
	  eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 > 0){
	  eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::vec Lambda01 = getLambda0WB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0WB(y1, a2, k2);
	arma::vec Lambda03 = getLambda0WB(y2-y1, a3, k3);
	arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

	arma::vec logdiff = arma::log(y2-y1);
	logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.

	arma::vec temp_scorevec(p1+p2+p3+7,arma::fill::zeros);
	//k1 (what ina calls u5)
	temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	//a1 (what ina calls u6)
	temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) -
									exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//k2 (what ina calls u7)
	temp_scorevec(2) = arma::accu( (1-delta1) % delta2 -
									commonVec % Lambda02 % arma::exp(h+eta2));
	//a2 (what ina calls u8)
	temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) -
									commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );
	//k3 (what ina calls u9)
	temp_scorevec(4) = arma::accu( delta1 % delta2 -
									commonVec % Lambda03 % arma::exp(h+eta3));
	//a3 (what ina calls u10)
	temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * logdiff) -
									commonVec % Lambda03 % arma::exp(h+a3+eta3) % logdiff );
	//h (what ina calls u1)
	temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec));
	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	}
	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	}
	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
	}
  return -temp_scorevec;
}

//this is the gradient with the markov assumption

// [[Rcpp::export]]
arma::vec ngradWB_ID_frail_M(const arma::vec& para,
							 const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
							 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
	  eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 > 0){
	  eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 > 0){
	  eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::vec Lambda01 = getLambda0WB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0WB(y1, a2, k2);
	arma::vec Lambda03y2 = getLambda0WB(y2, a3, k3);
	arma::vec Lambda03y1 = getLambda0WB(y1, a3, k3);
	arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + (Lambda03y2 - Lambda03y1) % arma::exp(eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

	arma::vec temp_scorevec(p1+p2+p3+7,arma::fill::zeros);
	//k1 (what ina calls u5)
	temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	//a1 (what ina calls u6)
	temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) -
									exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//k2 (what ina calls u7)
	temp_scorevec(2) = arma::accu( (1-delta1) % delta2 -
									commonVec % Lambda02 % arma::exp(h+eta2));
	//a2 (what ina calls u8)
	temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) -
									commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );
	//k3 (what ina calls u9)
	temp_scorevec(4) = arma::accu( delta1 % delta2 -
									commonVec % (Lambda03y2 - Lambda03y1) % arma::exp(h+eta3));
	//a3 (what ina calls u10)
	temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) -
									commonVec % arma::exp(h+a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) );
	//h (what ina calls u1)
	temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec));
	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	}
	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	}
	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % (Lambda03y2 - Lambda03y1)  % arma::exp(h + eta3));
	}
  return -temp_scorevec;
}

*/

// [[Rcpp::export]]
arma::mat ngradWB_ID_frail_mat_SM(const arma::vec& para,
                                  const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
                                  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
  //define constants
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

  double k1 = para(0);
  double a1 = para(1);
  double k2 = para(2);
  double a2 = para(3);
  double k3 = para(4);
  double a3 = para(5);
  double h = para(6);

  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
  }

  arma::vec Lambda01 = getLambda0WB(y1, a1, k1);
  arma::vec Lambda02 = getLambda0WB(y1, a2, k2);
  arma::vec Lambda03 = getLambda0WB(y2-y1, a3, k3);
  arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);
  arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

  arma::vec logdiff = arma::log(y2-y1);
  logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.

  arma::mat temp_scoremat(n,p1+p2+p3+7,arma::fill::zeros);
  //k1 (what ina calls u5)
  temp_scoremat(arma::span(0,n-1),0) = delta1 - commonVec % Lambda01 % arma::exp(h + eta1);
  //a1 (what ina calls u6)
  temp_scoremat(arma::span(0,n-1),1) = delta1 % (1 + exp(a1) * arma::log(y1)) -
    exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1);
  //k2 (what ina calls u7)
  temp_scoremat(arma::span(0,n-1),2) = (1-delta1) % delta2 -
    commonVec % Lambda02 % arma::exp(h+eta2);
  //a2 (what ina calls u8)
  temp_scoremat(arma::span(0,n-1),3) = (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) -
    commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1);
  //k3 (what ina calls u9)
  temp_scoremat(arma::span(0,n-1),4) = delta1 % delta2 -
    commonVec % Lambda03 % arma::exp(h+eta3);
  //a3 (what ina calls u10)
  temp_scoremat(arma::span(0,n-1),5) = delta1 % delta2 % (1 + exp(a3) * logdiff) -
    commonVec % Lambda03 % arma::exp(h+a3+eta3) % logdiff;
  //h (what ina calls u1)
  temp_scoremat(arma::span(0,n-1),6) = exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec);
  //beta1 (what ina calls u2)
  if(p1 > 0){
    temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)) = X1;
    temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)).each_col() %= (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
  }
  //beta2 (what ina calls u3)
  if(p2 > 0){
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2;
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)).each_col() %= ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
  }
  //beta3 (what ina calls u4)
  if(p3 > 0){
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1))  = X3;
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)).each_col() %= (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
  }
  return -temp_scoremat;
}


// [[Rcpp::export]]
arma::mat ngradWB_ID_frail_mat_M(const arma::vec& para,
                                 const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
                                 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
  //define constants
  int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

  double k1 = para(0);
  double a1 = para(1);
  double k2 = para(2);
  double a2 = para(3);
  double k3 = para(4);
  double a3 = para(5);
  double h = para(6);

  //define linear predictors
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n,arma::fill::zeros);
  if(p1 > 0){
    eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
  }
  if(p2 > 0){
    eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
  }
  if(p3 > 0){
    eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
  }

  arma::vec Lambda01 = getLambda0WB(y1, a1, k1);
  arma::vec Lambda02 = getLambda0WB(y1, a2, k2);
  arma::vec Lambda03y2 = getLambda0WB(y2, a3, k3);
  arma::vec Lambda03y1 = getLambda0WB(y1, a3, k3);
  arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + (Lambda03y2 - Lambda03y1) % arma::exp(eta3);
  arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

  arma::mat temp_scoremat(n,p1+p2+p3+7,arma::fill::zeros);
  //k1 (what ina calls u5)
  temp_scoremat(arma::span(0,n-1),0) = delta1 - commonVec % Lambda01 % arma::exp(h + eta1);
  //a1 (what ina calls u6)
  temp_scoremat(arma::span(0,n-1),1) = delta1 % (1 + exp(a1) * arma::log(y1)) -
    exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1);
  //k2 (what ina calls u7)
  temp_scoremat(arma::span(0,n-1),2) = (1-delta1) % delta2 -
    commonVec % Lambda02 % arma::exp(h+eta2);
  //a2 (what ina calls u8)
  temp_scoremat(arma::span(0,n-1),3) = (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) -
    commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1);
  //k3 (what ina calls u9)
  temp_scoremat(arma::span(0,n-1),4) = delta1 % delta2 -
    commonVec % (Lambda03y2 - Lambda03y1) % arma::exp(h+eta3);
  //a3 (what ina calls u10)
  temp_scoremat(arma::span(0,n-1),5) = delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) -
    commonVec % arma::exp(h+a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) ;
  //h (what ina calls u1)
  temp_scoremat(arma::span(0,n-1),6) = exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec);
  //beta1 (what ina calls u2)
  if(p1 > 0){
    temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)) = X1;
    temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)).each_col() %= (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
  }
  //beta2 (what ina calls u3)
  if(p2 > 0){
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2;
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)).each_col() %= ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
  }
  //beta3 (what ina calls u4)
  if(p3 > 0){
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1))  = X3;
    temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)).each_col() %= (delta1 % delta2 - commonVec % (Lambda03y2 - Lambda03y1) % arma::exp(h + eta3));
  }

  return -temp_scoremat;
}


//this is the hessian with the semi-markov assumption

// [[Rcpp::export]]
arma::mat nhessWB_ID_frail_SM(const arma::vec& para,
							  const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
							  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
	  eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 > 0){
	  eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 > 0){
	  eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::vec Lambda01 = getLambda0WB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0WB(y1, a2, k2);
	arma::vec Lambda03 = getLambda0WB(y2-y1, a3, k3);
	arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03 % arma::exp(eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

	arma::vec logdiff = arma::log(y2-y1);
	logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.

	arma::vec dAVecdk1 = Lambda01 % arma::exp(eta1);
	arma::vec dAVecdk2 = Lambda02 % arma::exp(eta2);
	arma::vec dAVecdk3 = Lambda03 % arma::exp(eta3);

	arma::vec dAVecda1 = Lambda01 % arma::exp(a1 + eta1) % arma::log(y1);
	arma::vec dAVecda2 = Lambda02 % arma::exp(a2 + eta2) % arma::log(y1);
	arma::vec dAVecda3 = Lambda03 % arma::exp(a3 + eta3) % logdiff;

	arma::vec dcommonVecdk1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk1;
	arma::vec dcommonVecdk2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk2;
	arma::vec dcommonVecdk3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk3;

	arma::vec dcommonVecda1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda1;
	arma::vec dcommonVecda2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda2;
	arma::vec dcommonVecda3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda3;

	arma::mat temp_hessmat(p1+p2+p3+7,p1+p2+p3+7,arma::fill::zeros);

	//======k1 section=======//
	//temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));

	//k1k1
	temp_hessmat(0,0) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	//k1a1
	temp_hessmat(0,1) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	temp_hessmat(1,0) = temp_hessmat(0,1);
	//k1k2
	temp_hessmat(0,2) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	temp_hessmat(2,0) = temp_hessmat(0,2);
	//k1a2
	temp_hessmat(0,3) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	temp_hessmat(3,0) = temp_hessmat(0,3);
	//k1k3
	temp_hessmat(0,4) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	temp_hessmat(4,0) = temp_hessmat(0,4);
	//k1a3
	temp_hessmat(0,5) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	temp_hessmat(5,0) = temp_hessmat(0,5);
	//k1h
	temp_hessmat(6,0) = arma::accu( dAVecdk1 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecdk1 + commonVec % dAVecdk1 ));
	temp_hessmat(0,6) = temp_hessmat(6,0);

	//==========k2 section==========//
	//temp_scorevec(2) = arma::accu( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h+eta2));

	//k2k2
	temp_hessmat(2,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02 ) );
	//k2a1
	temp_hessmat(1,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	temp_hessmat(2,1) = temp_hessmat(1,2);
	//k2a2
	temp_hessmat(3,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	temp_hessmat(2,3) = temp_hessmat(3,2);
	//k2k3
	temp_hessmat(4,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	temp_hessmat(2,4) = temp_hessmat(4,2);
	//k2a3
	temp_hessmat(5,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	temp_hessmat(2,5) = temp_hessmat(5,2);

	//===========k3 section==========//
	//temp_scorevec(4) = arma::accu( delta1 % delta2 - commonVec % Lambda03 % arma::exp(h+eta3));

	//k3k3
	temp_hessmat(4,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03 + commonVec % Lambda03 ) );
	//k3a1
	temp_hessmat(1,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03 ) );
	temp_hessmat(4,1) = temp_hessmat(1,4);
	//k3a2
	temp_hessmat(3,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03 ) );
	temp_hessmat(4,3) = temp_hessmat(3,4);
	//k3a3
	temp_hessmat(5,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03 + exp(a3) * commonVec % Lambda03 % logdiff ) );
	temp_hessmat(4,5) = temp_hessmat(5,4);

	//===========a1 section==========//
	//temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) - exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//rewritten as 		 arma::accu( delta1 + exp(a1) * delta1 % arma::log(y1) - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * commonVec % Lambda01) );

	//a1a1
	temp_hessmat(1,1) = arma::accu( exp(a1) * delta1 % arma::log(y1) -
									arma::exp(h+eta1) % arma::log(y1) % ( exp(a1) * commonVec % Lambda01 +
																		  exp(a1) * dcommonVecda1 % Lambda01 +
																		  exp(2 * a1) * commonVec % Lambda01 % arma::log(y1) ) );
	//a1a2
	temp_hessmat(1,3) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda2 % Lambda01) );
	temp_hessmat(3,1) = temp_hessmat(1,3);
	//a1a3
	temp_hessmat(1,5) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda3 % Lambda01) );
	temp_hessmat(5,1) = temp_hessmat(1,5);

	//===========a2 section==========//
	//temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );

	//rewritten as arma::accu( (1-delta1) % delta2 + exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02) );
	//a2a2
	temp_hessmat(3,3) = arma::accu( exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02 + exp(a2) * dcommonVecda2 % Lambda02 + exp(2 * a2) * commonVec % Lambda02 % arma::log(y1)) );
	//a2a3
	temp_hessmat(3,5) = arma::accu( - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * dcommonVecda3 % Lambda02) );
	temp_hessmat(5,3) = temp_hessmat(3,5);

	//=========a3 section==========//
	//temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * logdiff) - commonVec % Lambda03 % arma::exp(h+a3+eta3) % logdiff );

	//a3a3
	temp_hessmat(5,5) = arma::accu( exp(a3) * delta1 % delta2 % logdiff - arma::exp(h+eta3) % logdiff % (exp(a3) * commonVec % Lambda03 + exp(a3) * dcommonVecda3 % Lambda03 + exp(2 * a3) * commonVec % Lambda03 % logdiff) );

	//==========h section==========//
	//temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2 / (1+exp(h)) + arma::log(1+exp(h) * AVec)/exp(2*h) - commonVec % AVec) );

	//hh
	// from wolfram alpha (as an alternate derivation, with delta1=a and delta2=b):
	// https://www.wolframalpha.com/input/?i=derivative+of+e%5Eh+*+(+a*b+%2F+(1%2Be%5Eh)+%2B+log(1+%2B+e%5Eh+*+A)+%2F+(e%5E(2*h))+-+(e%5E(-h)+%2B+a+%2B+b)+*+A+%2F+(1+%2B+e%5Eh+*+A))+with+respect+to+h
	temp_hessmat(6,6) = arma::accu( (-delta1 + 2 * AVec - delta2) / (AVec * exp(h)+1)
								  + (delta1 - AVec + delta2) / (arma::pow(AVec * exp(h) + 1 , 2))
								  + delta1 % delta2 / (exp(h) + 1)
								  - delta1 % delta2 / (pow(exp(h) + 1, 2))
								  - exp(-h) * arma::log1p(AVec * exp(h))
								  );
	//ha1
	temp_hessmat(6,1) = arma::accu( dAVecda1 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecda1 + commonVec % dAVecda1 ));
	temp_hessmat(1,6) = temp_hessmat(6,1);
	//hk2
	temp_hessmat(6,2) = arma::accu( dAVecdk2 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecdk2 + commonVec % dAVecdk2 ));
	temp_hessmat(2,6) = temp_hessmat(6,2);
	//ha2
	temp_hessmat(6,3) = arma::accu( dAVecda2 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecda2 + commonVec % dAVecda2 ));
	temp_hessmat(3,6) = temp_hessmat(6,3);
	//hk3
	temp_hessmat(6,4) = arma::accu( dAVecdk3 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecdk3 + commonVec % dAVecdk3 ));
	temp_hessmat(4,6) = temp_hessmat(6,4);
	//ha3
	temp_hessmat(6,5) = arma::accu( dAVecda3 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecda3 + commonVec % dAVecda3 ));
	temp_hessmat(5,6) = temp_hessmat(6,5);

	//==========beta1 section=======//
	//temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	if(p1 > 0){
		//beta1k1
		temp_hessmat(arma::span(7,7 + p1 - 1),0) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	    temp_hessmat(0,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),0).t();
		//beta1a1
		temp_hessmat(arma::span(7,7 + p1 - 1),1) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	    temp_hessmat(1,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),1).t();
		//beta1k2
		temp_hessmat(arma::span(7,7 + p1 - 1),2) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	    temp_hessmat(2,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),2).t();
		//beta1a2
		temp_hessmat(arma::span(7,7 + p1 - 1),3) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	    temp_hessmat(3,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),3).t();
		//beta1k3
		temp_hessmat(arma::span(7,7 + p1 - 1),4) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	    temp_hessmat(4,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),4).t();
		//beta1a3
		temp_hessmat(arma::span(7,7 + p1 - 1),5) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	    temp_hessmat(5,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),5).t();
	    //beta1h
	    temp_hessmat(arma::span(7,7 + p1 - 1),6) = X1.t() * ( - Lambda01 % arma::exp(h + eta1) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),6).t();

	    //beta1beta1
	    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7,7 + p1 - 1)) = X1.t() * (X1.each_col() % ( - (Lambda01 % arma::exp(h + eta1) % ( commonVec + dcommonVecdk1 )))); //computes sum of w_i * x_ix_i^T

		if(p2 > 0){
		    //beta1beta2
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X1.t() * (X2.each_col() % ( - dcommonVecdk2 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)).t();
		}

		if(p3 > 0){
		    //beta1beta3
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X1.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//===========beta2 section=========//
    //temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	if(p2 > 0){
		//beta2k1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk1 % Lambda02 ) );
	    temp_hessmat(0,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0).t();
		//beta2a1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	    temp_hessmat(1,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1).t();
		//beta2k2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02) );
	    temp_hessmat(2,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2).t();
		//beta2a2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	    temp_hessmat(3,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3).t();
		//beta2k3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	    temp_hessmat(4,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4).t();
		//beta2a3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	    temp_hessmat(5,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5).t();
	    //beta2h
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6) = X2.t() * ( - Lambda02 % arma::exp(h + eta2) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6).t();

	    //beta2beta2
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X2.t() * (X2.each_col() % ( - (Lambda02 % arma::exp(h + eta2) % ( commonVec + dcommonVecdk2 )))); //computes sum of w_i * x_ix_i^T
	    if(p3 > 0){
		    //beta2beta3
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X2.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda02 % arma::exp(h + eta2))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//==========beta3 section=========//
    //temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
	if(p3 > 0){
		//beta3k1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk1 % Lambda03 ) );
	    temp_hessmat(0,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0).t();
		//beta3a1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03 ) );
	    temp_hessmat(1,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1).t();
		//beta3k2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk2 % Lambda03 ) );
	    temp_hessmat(2,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2).t();
		//beta3a2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03 ) );
	    temp_hessmat(3,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3).t();
		//beta3k3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03  + commonVec % Lambda03 ) );
	    temp_hessmat(4,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4).t();
		//beta3a3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03 + exp(a3) * commonVec % Lambda03 % logdiff ) );
	    temp_hessmat(5,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5).t();
	    //beta3h
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6) = X3.t() * ( - Lambda03 % arma::exp(h + eta3) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6).t();

	    //beta3beta3
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X3.t() * (X3.each_col() % ( - (Lambda03 % arma::exp(h + eta3) % ( commonVec + dcommonVecdk3 )))); //computes sum of w_i * x_ix_i^T
	}

	return -temp_hessmat;
}








//this is the hessian with the markov assumption

// [[Rcpp::export]]
arma::mat nhessWB_ID_frail_M(const arma::vec& para,
							 const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
							 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){

	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1(n,arma::fill::zeros);
	arma::vec eta2(n,arma::fill::zeros);
	arma::vec eta3(n,arma::fill::zeros);
	if(p1 > 0){
	  eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 > 0){
	  eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 > 0){
	  eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::vec Lambda01 = getLambda0WB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0WB(y1, a2, k2);
	arma::vec Lambda03y2 = getLambda0WB(y2, a3, k3);
	arma::vec Lambda03y1 = getLambda0WB(y1, a3, k3);
	arma::vec Lambda03diff = Lambda03y2 - Lambda03y1;
	arma::vec AVec = Lambda01 % arma::exp(eta1) + Lambda02 % arma::exp(eta2) + Lambda03diff % arma::exp(eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

	arma::vec dAVecdk1 = Lambda01 % arma::exp(eta1);
	arma::vec dAVecdk2 = Lambda02 % arma::exp(eta2);
	arma::vec dAVecdk3 = Lambda03diff % arma::exp(eta3);

	arma::vec dAVecda1 = Lambda01 % arma::exp(a1 + eta1) % arma::log(y1);
	arma::vec dAVecda2 = Lambda02 % arma::exp(a2 + eta2) % arma::log(y1);
	arma::vec dAVecda3 = arma::exp(a3 + eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1));

	arma::vec dcommonVecdk1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk1;
	arma::vec dcommonVecdk2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk2;
	arma::vec dcommonVecdk3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk3;

	arma::vec dcommonVecda1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda1;
	arma::vec dcommonVecda2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda2;
	arma::vec dcommonVecda3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda3;

	arma::mat temp_hessmat(p1+p2+p3+7,p1+p2+p3+7,arma::fill::zeros);

	//===========k1 section==========//
	//temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));

	//k1k1
	temp_hessmat(0,0) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	//k1a1
	temp_hessmat(0,1) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	temp_hessmat(1,0) = temp_hessmat(0,1);
	//k1k2
	temp_hessmat(0,2) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	temp_hessmat(2,0) = temp_hessmat(0,2);
	//k1a2
	temp_hessmat(0,3) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	temp_hessmat(3,0) = temp_hessmat(0,3);
	//k1k3
	temp_hessmat(0,4) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	temp_hessmat(4,0) = temp_hessmat(0,4);
	//k1a3
	temp_hessmat(0,5) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	temp_hessmat(5,0) = temp_hessmat(0,5);
	//k1h
	temp_hessmat(6,0) = arma::accu( dAVecdk1 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecdk1 + commonVec % dAVecdk1 ));
	temp_hessmat(0,6) = temp_hessmat(6,0);

	//===========k2 section===========//
	//temp_scorevec(2) = arma::accu( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h+eta2));

	//k2k2
	temp_hessmat(2,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02 ) );
	//k2a1
	temp_hessmat(1,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	temp_hessmat(2,1) = temp_hessmat(1,2);
	//k2a2
	temp_hessmat(3,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	temp_hessmat(2,3) = temp_hessmat(3,2);
	//k2k3
	temp_hessmat(4,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	temp_hessmat(2,4) = temp_hessmat(4,2);
	//k2a3
	temp_hessmat(5,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	temp_hessmat(2,5) = temp_hessmat(5,2);

	//===========k3 section===========//
	//temp_scorevec(4) = arma::accu( delta1 % delta2 - commonVec % Lambda03 % arma::exp(h+eta3));

	//k3k3
	temp_hessmat(4,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03diff + commonVec % Lambda03diff ) );
	//k3a1
	temp_hessmat(1,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03diff ) );
	temp_hessmat(4,1) = temp_hessmat(1,4);
	//k3a2
	temp_hessmat(3,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03diff ) );
	temp_hessmat(4,3) = temp_hessmat(3,4);
	//k3a3
	temp_hessmat(5,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03diff + exp(a3) * commonVec % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) ) );
	temp_hessmat(4,5) = temp_hessmat(5,4);

	//==========a1 section===========//
	//temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) - exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//rewritten as 		 arma::accu( delta1 + exp(a1) * delta1 % arma::log(y1) - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * commonVec % Lambda01) );

	//a1a1
	temp_hessmat(1,1) = arma::accu( exp(a1) * delta1 % arma::log(y1) -
									arma::exp(h+eta1) % arma::log(y1) % ( exp(a1) * commonVec % Lambda01 +
																		  exp(a1) * dcommonVecda1 % Lambda01 +
																		  exp(2 * a1) * commonVec % Lambda01 % arma::log(y1) ) );
	//a1a2
	temp_hessmat(1,3) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda2 % Lambda01) );
	temp_hessmat(3,1) = temp_hessmat(1,3);
	//a1a3
	temp_hessmat(1,5) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda3 % Lambda01) );
	temp_hessmat(5,1) = temp_hessmat(1,5);

	//==========a2 section============//
	//temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );

	//rewritten as arma::accu( (1-delta1) % delta2 + exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02) );
	//a2a2
	temp_hessmat(3,3) = arma::accu( exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02 + exp(a2) * dcommonVecda2 % Lambda02 + exp(2 * a2) * commonVec % Lambda02 % arma::log(y1)) );
	//a2a3
	temp_hessmat(3,5) = arma::accu( - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * dcommonVecda3 % Lambda02) );
	temp_hessmat(5,3) = temp_hessmat(3,5);

	//=========a3 section========//
	//temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) -
	//								commonVec % arma::exp(h+a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) );

	//a3a3
	temp_hessmat(5,5) = arma::accu( exp(a3) * delta1 % delta2 % arma::log(y2)
									- arma::exp(h+eta3) % ( exp(a3) * commonVec % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1))
														    + exp(a3) * dcommonVecda3 % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1))
														    + exp(2 * a3) * commonVec % (Lambda03y2 % arma::log(y2) % arma::log(y2) - Lambda03y1 % arma::log(y1) % arma::log(y1)) ) ) ;

	//==========h section===========//
	//temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2 / (1+exp(h)) + arma::log(1+exp(h) * AVec)/exp(2*h) - commonVec % AVec) );

	//hh
	// from wolfram alpha (as an alternate derivation, with delta1=a and delta2=b):
	// https://www.wolframalpha.com/input/?i=derivative+of+e%5Eh+*+(+a*b+%2F+(1%2Be%5Eh)+%2B+log(1+%2B+e%5Eh+*+A)+%2F+(e%5E(2*h))+-+(e%5E(-h)+%2B+a+%2B+b)+*+A+%2F+(1+%2B+e%5Eh+*+A))+with+respect+to+h
	temp_hessmat(6,6) = arma::accu( (-delta1 + 2 * AVec - delta2) / (AVec * exp(h)+1)
								  + (delta1 - AVec + delta2) / (arma::pow(AVec * exp(h) + 1 , 2))
								  + delta1 % delta2 / (exp(h) + 1)
								  - delta1 % delta2 / (pow(exp(h) + 1, 2))
								  - exp(-h) * arma::log1p(AVec * exp(h))
								  );
	//ha1
	temp_hessmat(6,1) = arma::accu( dAVecda1 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecda1 + commonVec % dAVecda1 ));
	temp_hessmat(1,6) = temp_hessmat(6,1);
	//hk2
	temp_hessmat(6,2) = arma::accu( dAVecdk2 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecdk2 + commonVec % dAVecdk2 ));
	temp_hessmat(2,6) = temp_hessmat(6,2);
	//ha2
	temp_hessmat(6,3) = arma::accu( dAVecda2 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecda2 + commonVec % dAVecda2 ));
	temp_hessmat(3,6) = temp_hessmat(6,3);
	//hk3
	temp_hessmat(6,4) = arma::accu( dAVecdk3 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecdk3 + commonVec % dAVecdk3 ));
	temp_hessmat(4,6) = temp_hessmat(6,4);
	//ha3
	temp_hessmat(6,5) = arma::accu( dAVecda3 / (1 + exp(h)*AVec) -
							exp(h) * (AVec % dcommonVecda3 + commonVec % dAVecda3 ));
	temp_hessmat(5,6) = temp_hessmat(6,5);

	//=========beta1 section===========//
	//temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	if(p1 > 0){
		//beta1k1
		temp_hessmat(arma::span(7,7 + p1 - 1),0) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	    temp_hessmat(0,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),0).t();
		//beta1a1
		temp_hessmat(arma::span(7,7 + p1 - 1),1) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	    temp_hessmat(1,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),1).t();
		//beta1k2
		temp_hessmat(arma::span(7,7 + p1 - 1),2) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	    temp_hessmat(2,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),2).t();
		//beta1a2
		temp_hessmat(arma::span(7,7 + p1 - 1),3) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	    temp_hessmat(3,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),3).t();
		//beta1k3
		temp_hessmat(arma::span(7,7 + p1 - 1),4) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	    temp_hessmat(4,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),4).t();
		//beta1a3
		temp_hessmat(arma::span(7,7 + p1 - 1),5) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	    temp_hessmat(5,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),5).t();
	    //beta1h
	    temp_hessmat(arma::span(7,7 + p1 - 1),6) = X1.t() * ( - Lambda01 % arma::exp(h + eta1) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),6).t();

	    //beta1beta1
	    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7,7 + p1 - 1)) = X1.t() * (X1.each_col() % ( - (Lambda01 % arma::exp(h + eta1) % ( commonVec + dcommonVecdk1 )))); //computes sum of w_i * x_ix_i^T

	    if(p2 > 0){
		    //beta1beta2
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X1.t() * (X2.each_col() % ( - dcommonVecdk2 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)).t();
		}

	    if(p3 > 0){
		//beta1beta3
	    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
	    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X1.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//=========beta2 section==========//
    //temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	if(p2 > 0){
		//beta2k1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk1 % Lambda02 ) );
	    temp_hessmat(0,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0).t();
		//beta2a1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	    temp_hessmat(1,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1).t();
		//beta2k2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02) );
	    temp_hessmat(2,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2).t();
		//beta2a2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	    temp_hessmat(3,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3).t();
		//beta2k3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	    temp_hessmat(4,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4).t();
		//beta2a3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	    temp_hessmat(5,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5).t();
	    //beta2h
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6) = X2.t() * ( - Lambda02 % arma::exp(h + eta2) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6).t();

	    //beta2beta2
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X2.t() * (X2.each_col() % ( - (Lambda02 % arma::exp(h + eta2) % ( commonVec + dcommonVecdk2 )))); //computes sum of w_i * x_ix_i^T

		if(p3 > 0){
		    //beta2beta3
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X2.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda02 % arma::exp(h + eta2))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//=========beta3 section===========//
	//temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % (Lambda03y2 - Lambda03y1)  % arma::exp(h + eta3));
    //temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
	if(p3 > 0){
		//beta3k1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk1 % Lambda03diff ) );
	    temp_hessmat(0,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0).t();
		//beta3a1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03diff ) );
	    temp_hessmat(1,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1).t();
		//beta3k2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk2 % Lambda03diff ) );
	    temp_hessmat(2,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2).t();
		//beta3a2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03diff ) );
	    temp_hessmat(3,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3).t();
		//beta3k3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03diff  + commonVec % Lambda03diff ) );
	    temp_hessmat(4,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4).t();
		//beta3a3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03diff + exp(a3) * commonVec % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) ) );
	    temp_hessmat(5,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5).t();
	    //beta3h
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6) = X3.t() * ( - Lambda03diff % arma::exp(h + eta3) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6).t();

	    //beta3beta3
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X3.t() * (X3.each_col() % ( - (Lambda03diff % arma::exp(h + eta3) % ( commonVec + dcommonVecdk3 )))); //computes sum of w_i * x_ix_i^T
	}
	return -temp_hessmat;
}



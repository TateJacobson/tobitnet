#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

double normalCDF(double x){
    return (std::erfc(-x/std::sqrt(2))/2);
}

double normalPDF(double x){
	static const double inv_sqrt_2pi = 0.3989422804014327;
	return(inv_sqrt_2pi * std::exp(-0.5 * x * x));
}

double g(double x){
    double gout;
	if(x < -37){
		gout = -x;
	}else{
		gout = normalPDF(x)/normalCDF(x);
	}
    return(gout);
}

//[[Rcpp::export]]
NumericVector LprimeC(NumericVector xj, NumericVector y, LogicalVector d, NumericVector r, double gamma){
	int n = y.size();
	NumericVector out(n);
	for(int i = 0; i < n; i++){
		if(d[i]){
			out[i] = -(gamma*y[i] - r[i])*xj[i];
		}else{
			out[i] = g(-r[i])*xj[i];
		}
	}
	return out;
}

//[[Rcpp::export]]
NumericVector logL1(NumericVector y, LogicalVector d, NumericVector r, double gamma){
	int n = y.size();
	NumericVector out(n);
	for(int i = 0; i < n; i++){
		if(d[i]){
			out[i] = std::log(gamma) -(gamma*y[i] - r[i])*(gamma*y[i] - r[i])/2;
		}else{
			out[i] = std::log(normalCDF(-r[i])) ;
		}
	}
	return out;
}

double soft(double z, double t){
	double sgn = (z > 0) - (z < 0);
	double out = std::max(0.0, std::abs(z) - t)*sgn;
	return out;
}

bool KKT_check(double deltaj, NumericVector xj, NumericVector y, LogicalVector d, NumericVector r, double gamma, double l1, double l2, double tol = 1e-3){
	if(deltaj > 0){
        double check = mean( LprimeC(xj, y, d, r, gamma) ) + l2*deltaj + l1;
        if( -tol < check && check < tol ) return(true);	
    } else if (deltaj < 0){
        double check = mean( LprimeC(xj, y, d, r, gamma) ) + l2*deltaj - l1;
        if( -tol < check && check < tol ) return(true);		
    } else if (deltaj == 0){
        double check = - mean( LprimeC(xj, y, d, r, gamma) ); 
        if( - l1 - tol  < check && check < l1 + tol ) return(true); 
    }
    return(false);
}

NumericVector matProd(NumericMatrix x, NumericVector b){
	int n = x.nrow();
	NumericVector out(n);
	
	for(int i = 0; i<n; i++){
		out(i) = sum(x(i,_)*b);
	}
	
	return out;
}

List standardizeC(NumericMatrix x){
	int p = x.ncol();
	NumericVector colmeans(p);
	NumericVector colsds(p);
	
	for(int j = 0; j < p; j++){
		colmeans(j) = mean( x(_,j) );
		x(_,j) = x(_,j) - colmeans(j);
		colsds(j) = sd( x(_,j) );
		x(_,j) = x(_,j)/colsds(j);
	}
	
	return List::create(
	Named("colmeans") = colmeans,
	Named("colsds") = colsds
	);
}

//[[Rcpp::export]]
List tobitnet_innerC(NumericMatrix xin, NumericVector yin, NumericVector cin, double lambda1, double lambda2, NumericVector pf1, NumericVector pf2, NumericVector delta_init, double delta_0_init = 0, double gamma_init = 1, double eps = 1e-7, bool standardize = true, double maxit = 1e6){
	NumericMatrix x = clone(xin);
	NumericVector y = clone(yin);
	NumericVector c = clone(cin);
	int n = y.size();
	int p = x.ncol();
	
	LogicalVector d(n);
	NumericVector r(n);
	LogicalVector active_set_bool = rep(true, p);
	LogicalVector current_active_set_bool = rep(true, p);
	
	NumericVector colmeans(p);
	NumericVector colsds(p);
	
	double delta_0_current = delta_0_init;
	NumericVector delta_current(p);
    NumericVector delta_step(p);
    NumericVector delta_new(p);
	
	if(is_finite(c)[0]){
		for(int i = 0; i < n; i++){
			y(i) = y(i) - c(0);
			d(i) = (y(i) > 0);
		}
    } else {
		d = rep(true, n);
    }
	
	if(p != 0){
        if(standardize){
			List xL = standardizeC(x);
			
			colmeans = xL["colmeans"];
			colsds = xL["colsds"];
        }
		
        //Initialized delta_0, delta, gamma (pass these in for the outer loop to implement warms starts)
        delta_current = delta_init;
        delta_step = rep(0.0, p);
        delta_new = rep(0.0, p);
        
		//updated often throughout loop
        r = matProd(x, delta_current) + rep(delta_0_current, n); 
		
    } else {
		standardize = false;
        r = rep(delta_0_current, n);
    }
	
	double ga = sum( pow(y, 2) );
	double ga2inv = 1/(2*ga);
	double gc = - sum(d);
    
    double gamma_current = gamma_init;
    double delta_0_step = 0.0;
    
	bool full_loop_check = false;
	
	NumericVector M(p);
	if(standardize){
		M = rep(1.0, p);
	}else{
		for(int j = 0; j < p; j++){
			M[j] = mean( pow(x(_,j), 2) );
		}
	}
	
	NumericVector updateDenom = 1/( M + lambda2*pf2 );
	NumericVector updatet = lambda1*pf1;
	
	for(int k_outer = 0 ; k_outer < maxit ; k_outer++){
        //Update delta_0
        delta_0_step = - mean( LprimeC(rep(1.0, n), y, d, r, gamma_current) );
        delta_0_current += delta_0_step;
        
        //Update r with new delta_0
        r = r + rep(delta_0_step, n);
        
        //Update delta
        for(int j = 0; j < p; j++){
			
			if(active_set_bool(j)){
				double z = M(j)*delta_current(j) - mean( LprimeC(x(_,j), y, d, r, gamma_current) );
				delta_new(j) = soft(z, updatet(j))*updateDenom(j); 
				
				delta_step(j) = delta_new(j) - delta_current(j);
				delta_current(j) = delta_new(j);
				
				//Update r with new delta_j
				r = r + delta_step(j)*x(_,j);	
				
				if( delta_current(j) == 0 && delta_step(j) == 0 ){				
					active_set_bool(j) = false;
				}
				
			}
			
        }
        
        //Update gamma
        double gb = - sum(y*r);
        gamma_current = (-gb + std::sqrt( pow(gb, 2) - 4*ga*gc))*ga2inv; 
        
        double max_delta2 = max( pow(delta_step, 2) );
        if( std::max(max_delta2, pow(delta_0_step, 2) ) < eps || full_loop_check ) {
			if( !full_loop_check ){
                current_active_set_bool = clone(active_set_bool);
				active_set_bool = rep(true, p); 
                full_loop_check = true;
            } else {
                if( is_true( all( active_set_bool == current_active_set_bool ) ) ){
					break;
                } else {
					full_loop_check = false;
				}
          
            }
        
		}
		
	}
	
	bool KKT = true;
	if(! KKT_check(delta_0_current, rep(1.0,n), y, d, r, gamma_current, 0, 0) ){
		KKT = false;
	}
	if(KKT){
		for(int j = 0; j < p; j++){
			if( ! KKT_check(delta_current(j), x(_,j), y, d, r, gamma_current, pf1(j)*lambda1, pf2(j)*lambda2) ){
				KKT = false;
				break;
			}	
		}
	}
	
	//Create beta and sigma, with proper rescaling if standardize == T (note that we do not rescale d0 and delta since we only use them for warm starts)
	NumericVector beta(p);
	
	double beta0 = delta_0_current/gamma_current;
	for(int j = 0; j < p; j++){
		beta(j) = delta_current(j)/gamma_current; 
	}
	
	if(standardize){
		for(int j = 0; j < p; j++){
				beta(j) = beta(j)/colsds(j);
			}
		beta0 -= sum(beta*colmeans);
	}
	
	if(is_finite(c)[0]){ beta0 += c(0); }
	
	double sigma = 1/gamma_current;
	
	return List::create(
	Named("d0") = delta_0_current,
	Named("delta") = delta_current,
	Named("gamma") = gamma_current,
	Named("b0") = beta0,
	Named("beta") = beta,
	Named("sigma") = sigma,
	Named("KKT") = KKT
	);
}

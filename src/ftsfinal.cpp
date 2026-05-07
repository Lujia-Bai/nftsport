#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace arma;
using namespace Rcpp;

#define TOLERANCE 1.0e-200
#define NA_VALUE R_NaReal
#define SQRT2 1.414214
//bootstrap

// [[Rcpp::export]]
double EpaK(double x)
  /*	Epanechnikov kernel	*/
{
  if(fabs(x) <= 1) return( 3*(1-x*x)/4 ); else return( 0 );
}

double EpaK0(double x1_0, double x2_0)
  /* v0	*/
{
  double x1,x2;
  if(x1_0 >=-1) x1 = x1_0;
  else x1 = -1;

  if(x2_0 <= 1) x2 = x2_0;
  else x2 = 1;

  return(0.75*((x2 - pow(x2,3)/3)-(x1 - pow(x1,3)/3)));
}

double EpaK1(double x1_0, double x2_0)
  /* v0	*/
{
  double x1,x2;
  if(x1_0 >=-1) x1 = x1_0;
  else x1 = -1;

  if(x2_0 <= 1) x2 = x2_0;
  else x2 = 1;

  return(0.75*((pow(x2,2)/2 - pow(x2,4)/4)-(pow(x1,2)/2 - pow(x1,4)/4)));
}

double EpaK2(double x1_0, double x2_0)
  /* v0	*/
{
  double x1,x2;
  if(x1_0 >=-1) x1 = x1_0;
  else x1 = -1;

  if(x2_0 <= 1) x2 = x2_0;
  else x2 = 1;

  return(0.75*((pow(x2,3)/3 - pow(x2,5)/5)-(pow(x1,3)/3 - pow(x1,5)/5)));
}



double db_EpaK(double x)
  /*	Epanechnikov kernel	*/
{
  return( 2*SQRT2*EpaK(SQRT2*x)- EpaK(x));
}



double psum(NumericVector e, double mean_e, int m, int j)
{
  double sjm;

  sjm = - m * mean_e;
  for(int k = 0; k <= m - 1; k++){
    sjm += e[j + k - 1];
    //j-1 becuase we take j as defined in the paper, starting from 1
  }
  return(sjm);
}

NumericVector all_psum(int m, int n, NumericVector e, bool mean = 1){
  double mean_e = 0.0;
  if(mean) mean_e = std::accumulate(e.begin(),e.end(),0.0)/n;
  NumericVector all_partial_sum(n - m + 1);

  for(int i = 0; i <  n - m + 1; i++){
    all_partial_sum[i] = psum(e, mean_e, m , i + 1); // all the S_{i,m}-m/n*S_n needed
  }
  return(all_partial_sum);
}


NumericVector Compute_kernel_vector(int n, double bw, int type){
  //type: 0 for naive kernel 1 for Jackknife 2 for Jackknife equivalent kernel
  NumericVector kernel_vector(n);

  if(type == 0){
    for(int i = 0; i < n * bw ;i++)
      kernel_vector[i] =  EpaK(((i + 0.)/n)/bw);//interval
  }else if(type == 1 || type == 2){
    for(int i = 0; i < n * bw ;i++)
      kernel_vector[i] =  db_EpaK(((i + 0.)/n)/bw);//interval
  }
  else{
    Rcout<<"Bad type!Please choose from 0,1,2"<<endl;
  }

  return(kernel_vector);
}
// [[Rcpp::export]]
arma::cube locLinSmootherC(double bw, arma::vec x, arma::cube y, bool db_kernel =0){
  //int n = y.size();
  int i, j;
  double s0, s1, s2, auxK, aux;
  double den;
  int n = y.n_rows;
  int N = y.n_cols;
  int p = y.n_slices;
  arma::cube t0(1, N, p), t1(1, N, p);
  arma::cube res(n, N, p);

  res.zeros();
  for(i = 0; i < n; i++){
    s0 = s1 = s2  = 0.;
    t0.zeros();
    t1.zeros();
    den = 0.;

    for(j = 0; j < n; j++){
      aux = (x[j]-x[i]) / bw;

      if(db_kernel)
      {//debias kernel
        auxK = db_EpaK(aux)/ bw;
      }else{
        auxK = EpaK(aux)/ bw;
      }

      t0 += y.row(j) * auxK;
      t1 += y.row(j) * auxK * aux;
      s0 += auxK;
      s1 += auxK * aux;
      s2 += auxK * aux * aux;
    }
    den = (s2*s0 - s1*s1);
    if(abs(den)> TOLERANCE){
      res.row(i) = y.row(i) - (s2*t0 - s1*t1)/den;
    }
    else
      res.row(i) = NA_VALUE;
  }
  return(res);
}




   //' @details
   //' The portmanteau test statistic for multivariate functional time series.\cr
   //' The time varying coefficients are estimated by
   //' \deqn{Q_n = \underset{1 \leq i \leq n, 1 \leq j \leq N}{\max}\sum_{k=1}^{s_n} ((n\tau) \hat {R}_{k}(i/n,j/N)- \phi^2)}
   //' where \eqn{\hat {R}_{k}(t, u)  =  |\hat {\bs \Gamma}_{k}(t, u)|^2/ | \hat {\bs \Gamma}_{0}(t, u)|^2}
   // [[Rcpp::export]]
   double ftsQ(double tau, int trim, arma::vec t, arma::cube e, int sn){

     int i, j, k;
     double auxK, aux, tmp;
     int n = e.n_rows;
     int N = e.n_cols;
     int p = e.n_slices;
     // mat t0(n, N);
     arma::mat tk(n, N);
     double Q;
     arma::mat temp(p,p), temp1(p,p);
     vec e1, e2;
     int m = n*tau + 1 ;
     int lb, ub;
     if(trim > sn){
       trim = sqrt(sn);
     }

     List result;


     tk.zeros();

     for(i = m + sn; i <= n-m; i++){
       for(int l = 0; l < N; l++){
         if(i-m-1 > sn)
           lb = i-m-1;
         else
           lb = sn;
         if(i+m+1 < n)
           ub = i+m+1;
         else ub = n;
         for(k = 1; k <= sn; k++ ){
           temp.zeros();
           for(j = lb; j < ub; j++){
             temp1.zeros();
             for(int q = sn; q < j - trim; q++){
               if(q >= j - sn && (q < n)){
                 aux = (t[q]-t[i]) / tau;
                 auxK = EpaK(aux);
                 e1 = vectorise(e.tube(q-k,l));
                 e2 = vectorise(e.tube(q,l));
                 temp1 += e1 * e2.t() * auxK;
               }
             }
             // if(j == lb){
             //  for(int q = sn; q < j + floor(n*a); q++){
             //      if(q >= j - floor(n*a) && (q < n)){
             //       aux = (t[q]-t[i]) / tau;
             //       auxK = EpaK(aux)/(2*n*a);
             //       e1 = vectorise(e.tube(q-k,l));
             //       e2 = vectorise(e.tube(q,l));
             //       temp1 += e1 * e2.t() * auxK;
             //     }
             //   }
             // }else{
             //     if(j + floor(n*a) - sn >= 0 && j+ floor(n*a) < n){
             //     aux = (t[j + floor(n*a)]- t[i]) / tau;
             //     auxK = EpaK(aux)/(2*n*a);
             //     e1 = vectorise(e.tube(j + floor(n*a) - k, l));
             //     e2 = vectorise(e.tube(j + floor(n*a), l));
             //     temp1 += e1 * e2.t() * auxK;
             //   }

             //   if(j - floor(n*a)-1-sn >= 0){
             //     aux = (t[j - floor(n*a) - 1]-t[i]) / tau;
             //     auxK = EpaK(aux)/(2*n*a);
             //     e1 = vectorise(e.tube(j - floor(n*a) - k - 1, l));
             //     e2 = vectorise(e.tube(j - floor(n*a) - 1, l));
             //     temp1 -= e1 * e2.t() * auxK;
             //   }
             // }
             // }

             aux = (t[j]-t[i]) / tau;
             auxK = EpaK(aux);
             e1 = vectorise(e.tube(j-k,l));
             e2 = vectorise(e.tube(j,l));
             temp += e1 * e2.t()* auxK * temp1.t() ;
           }

           tmp = trace(temp);
           // }

           tk(i, l) += tmp;

         }
       }
     }

     Q = max(abs(vectorise(tk)/sqrt(n*tau*(sn-trim))));

     return(Q);
   }


 //' @details
 //' The portmanteau test bootstrap multiplier for multivariate functional time series.\cr
 //' \deqn{\hat{\mathbf V}_i = (\hat {\mathbf V}_{i,\lfoor n\tau \rfloor } ,\cdots,\hat {\\mathbf V}_{i + n - 2\lfoor n\tau \rfloor , n-\lfoor n\tau \rfloor } )^{\top}}
 //' where \eqn{\hat {\mathbf V}_{i,l}:= K_{\tau}(i/n - l/n) (\mathrm{tr} (\hat{\mathbf U}_{i}(l/n, j/N))/|\hat{\boldsymbol \Gamma}_0(l/n, j/N)|^2, 1 \leq j \leq N)^{\top}}, and
 //' \eqn{\mathbf{ \mf U}_i(t, u) =\sum_{k=1}^{s_n}(\hat{\boldsymbol \phi}_{i,k}(u)\hat{\boldsymbol \psi}^{\top}_{i,k}(t,u) + \hat{\boldsymbol \psi}_{i,k}(t,u)\hat{\boldsymbol \phi}^{\top}_{i,k}(u) - \hat{\boldsymbol \phi}_{i,k}(u) \hat{\boldsymbol \phi}_{i,k}^{\top}(u)K_{\tau}(t_i -t))}
 // [[Rcpp::export]]
 arma::mat Qmultiplier(double tau, int trim, arma::vec t, arma::cube e, int sn){

   int n = e.n_rows;
   int N = e.n_cols;
   int p = e.n_slices;
   int lb, ub, m = 0;
   // mat U(n, N);
   List result;
   int k;

   if(trim > sn){
     trim = sqrt(sn);
   }

   arma::mat phi(p,p), psi(p,p), temp(p,p);
   // ,temp1(p,p), meanGk1(p,p), meanGk2(p,p);
   vec e1, e2;
   m = n*tau + 1;
   arma::mat V(2 * m, N  * ( n - 2 * m - sn+ 1));
   // ,V1(2 * m, N  * ( n - 2 * m - sn+ 1));
   double aux, auxi, auxK, auxKi;
   V.zeros();
   // U.zeros();
   // V1.zeros();



   for(int l = m + sn; l <= n - m; l++){
     for(int j = 0; j < N; j++){
       for(int i = 0; i < 2*m; i++){
         auxi = (i-m) / (n*tau);
         auxKi =  EpaK(auxi);
         // std::cout<< i <<endl;
         temp.zeros();
         // meanGk2.zeros();

         for(k = 1; k <= sn; k++ ){
           e1 = vectorise(e.tube(i + l - m - k, j));
           e2 = vectorise(e.tube(i + l - m, j));
           phi =  e1 * e2.t();
           psi.zeros();
           // if(l - m - 1 >  sn) lb = l - m - 1; else  lb = sn;
           if(i + l - m - sn - 1 >  sn) lb = i + l - m - sn - 1; else  lb = sn;
           if(i + l - m > n) ub = n; else  ub = i + l - m - trim;

           for(int h = lb; h < ub; h++){
             e1 = vectorise(e.tube(h - k, j));
             e2 = vectorise(e.tube(h, j));
             aux = (t[h] - t[l]) / tau;
             auxK = EpaK(aux);
             psi +=   e1 * e2.t() * auxK;
             // if(k == 1)
             //   meanGk2 += reshape(vectorise(Gk.tube(h,j)), p,p) * auxK;
           }

           // temp +=  phi*psi.t() + psi*phi.t()  - phi*phi.t()*auxKi;
           // temp +=  phi*psi.t() + psi*phi.t();
           temp +=  phi*psi.t();
         }
         // meanGk1 = reshape(vectorise(Gk.tube(i + l - m,j)), p,p);
         // temp1 = temp - meanGk1*psi.t() - phi*meanGk2.t() + meanGk1*meanGk2.t() - psi*meanGk1.t() - meanGk2*phi.t() + meanGk2*meanGk1.t()  + meanGk1*phi.t() + phi*meanGk1.t() - meanGk1*meanGk1.t();

         // if(std){
         //    U(l,j)  += auxKi * trace(temp)/G0(l,j);
         //    V(i, j + (l - m - sn) * N) += auxKi * trace(temp)/G0(l,j);
         //   //  V1(i, j + (l - m - sn) * N) += auxKi * trace(temp1)/G0(l,j);
         // }else{
         // U(l,j)  += auxKi * trace(temp) ;
         V(i, j + (l - m - sn) * N) += auxKi * trace(temp);
         // V1(i, j + (l - m - sn) * N) += auxKi * trace(temp1);
         // }


       }
     }
   }
   // result["V"] = V;
   // result["U"] = U;
   // result["V1"] = V1;
   return(V/sqrt(sn - trim)); //normalization
 }

 //' @details
 //' difference based multiplier
 // [[Rcpp::export]]
 arma::mat Qdist(arma::mat V, int n, double tau,  int N, int B, int sn,  int L = 5){

   arma::mat R(n, B);
   R.randn();

   int rowV = V.n_rows;
   int colV = V.n_cols;
   // std::cout << rowV-2*L+1 << endl;
   int rowS = rowV - 2*L+1;
   arma::mat S(rowS, colV);
   arma::mat output(1,B);
   arma::mat boot(n - rowV - sn + 1, B),  tmp(N, B);
   int col1, col2;
   double frac;

   S.zeros();
   boot.zeros();
   output.zeros();

   for(int i = sn; i < L+sn; i++){
     S.row(sn) += V.row(i+L) - V.row(i); //L:2*L-1 - 0:L-1
   }
   for(int j = 1 + sn; j < rowS; j++){
     S.row(j) =   S.row(j - 1) + V.row(j + 2 * L - 1) + V.row(j-1) - 2*V.row(j + L - 1);
   }
   S/= sqrt(2*L);
   // std::cout<<max(vectorise(S))<<endl;
   for(int l = 0; l <= n - rowV - sn; l++){
     tmp.zeros();
     for(int j = sn; j < rowS; j++){
       col1 = l * N;
       col2 = l * N + N - 1;
       tmp += S.submat(j, col1, j, col2).t() * R.row(l+j);
     }
     tmp = abs(tmp);
     boot.row(l) = max(tmp); // column max
   }
   //  std::cout<<max(vectorise(boot))<<endl;
   frac = sn* (rowS) /2;
   //  std::cout<<frac<<endl;
   output = max(boot)/sqrt(frac);
   return(output);
 }


////universal


 //' @details
 //' The extended portmanteau test bootstrap multiplier for multivariate functional time series.\cr
 //' \deqn{\hat{\mathbf V}_i = (\hat {\mathbf V}_{i,\lfoor n\tau \rfloor } ,\cdots,\hat {\\mathbf V}_{i + n - 2\lfoor n\tau \rfloor , n-\lfoor n\tau \rfloor } )^{\top}}
 //' where \eqn{\hat {\mathbf V}_{i,l}:= K_{\tau}(i/n - l/n) (\mathrm{tr} (\hat{\mathbf U}_{i}(l/n, j/N))/|\hat{\boldsymbol \Gamma}_0(l/n, j/N)|^2, 1 \leq j \leq N)^{\top}}, and
 //' \eqn{\mathbf{ \mf U}_i(t, u) =\sum_{k=1}^{s_n}(\hat{\boldsymbol \phi}_{i,k}(u)\hat{\boldsymbol \psi}^{\top}_{i,k}(t,u) + \hat{\boldsymbol \psi}_{i,k}(t,u)\hat{\boldsymbol \phi}^{\top}_{i,k}(u) - \hat{\boldsymbol \phi}_{i,k}(u) \hat{\boldsymbol \phi}_{i,k}^{\top}(u)K_{\tau}(t_i -t))}
 // [[Rcpp::export]]
 arma::mat Qmultiplierj_uni(double tau, int q, int trim, arma::vec t, arma::cube e, arma::cube dotmat, int sn){

   int n = e.n_rows;
   int N = e.n_cols;
   int p = e.n_slices;

   int lb, ub, m = 0;
   List result;
   int k;
   double temp = 0;

   if(trim > sn){
     trim = sqrt(sn);
   }

   arma::mat  psi(p,p);
   // vec e1, e2;
   // double e1, e2, phi, psi, temp;
   m = n*tau + 1;
   arma::mat V(2 * m, N  * ( n - 2 * m - sn+ 1));

   double aux, auxi, auxK, auxKi;
   V.zeros();


   for(int l = m + sn; l <= n - m; l++){
     for(int j = 0; j < N; j++){
       for(int i = 0; i < 2*m; i++){
         auxi = (i-m) / (n*tau);
         auxKi =  EpaK(auxi);
         // std::cout<< i <<endl;
         temp = 0;
         // meanGk2.zeros();

         for(k = 1; k <= sn; k++ ){
           psi.zeros();
           // if(l - m - 1 >  sn) lb = l - m - 1; else  lb = sn;
           if(i + l - m - sn - 1 >  sn) lb = i + l - m - sn - 1; else  lb = sn;
           if(i + l - m > n) ub = n; else  ub = i + l - m - trim;

           for(int h = lb; h < ub; h++){
             // e1 = vectorise(e.tube(h - k, j));
             // e2 = vectorise(e.tube(h, q));
             aux = (t[h] - t[l]) / tau;
             auxK = EpaK(aux);
             // Rcout<< "use dotmat" << i+l-m << i+l-m-h;
             temp += dotmat(i+l-m, i+l-m-h, q) * dotmat(i+l-m-k, i+l-m-h, j)* auxK;
             // psi +=   e1 * e2.t() * auxK;
           }
           // e1 = vectorise(e.tube(i + l - m - k, j));
           // e2 = vectorise(e.tube(i + l - m, q));
           // temp +=  arma::dot(e2, psi.t() * e1);
         }

         V(i, j + (l - m - sn) * N) += auxKi *  temp ;

       }
     }
   }

   return(V/sqrt(sn - trim)); //normalization
 }

// [[Rcpp::export]]
arma::cube compute_dotmat(arma::cube e, int sn){
  int n = e.n_rows;
  int N = e.n_cols;
  int p = e.n_slices;

 arma::cube dotmat(n, sn + 2, N);
  arma::vec ej_h(p), eq_h(p);


  for(int h = 0; h < N; h++){
    for (int j = sn; j < n; j++) {
      if(j == sn){
        for (int q = 0; q <= sn; q++) {
          ej_h = e.tube(j, h) ;
          eq_h = e.tube(j-q, h) ;
          dotmat(j, q, h) = dot(ej_h, eq_h);
        }
      }else{
        for (int q = 0; q <= sn + 1; q++) {
          ej_h = e.tube(j, h) ;
          eq_h = e.tube(j-q, h) ;
          dotmat(j, q, h) = dot(ej_h, eq_h);
        }
      }
    }
  }
  return dotmat;
}
// [[Rcpp::export]]
double ftsQcross_universal(double tau, int trim, arma::vec t, arma::cube e, arma::cube dotmat, int sn) {
  int n = e.n_rows;
  int N = e.n_cols;
  double temp = 0;

  arma::mat tk(n, N);
  double Q = 0.0;

  // Pre-calculate window parameters
  int m = n * tau + 1;
  if (trim > sn) trim = std::sqrt(sn);

  // Loop through each series 'h' (The "Cross" dimension)
  for (int h = 0; h < N; h++) {
    tk.zeros();

    for (int i = m + sn; i <= n - m; i++) {
      for (int l = 0; l < N; l++) {

        int lb = (i - m - 1 > sn) ? (i - m - 1) : sn;
        int ub = (i + m + 1 < n) ? (i + m + 1) : n;

        for (int k = 1; k <= sn; k++) {
          temp = 0;

          for (int j = lb; j < ub; j++) {
            double temp1 = 0;

            // ejk_l = e.tube(j - k, l);
            // ej_h = e.tube(j, h);
            // Internal summation over the lag window
            for (int q = sn; q < j - trim; q++) {
              if (q >= j - sn && q < n) {
                double auxK = EpaK((t[q] - t[i]) / tau);

                // Direct vector extraction from the cube
                // eqk_l = e.tube(q - k, l);
                // eq_h = e.tube(q, h);
                // temp1 +=  dot(ej_h , eq_h)  *  dot( ejk_l, eqk_l) * auxK;
                temp1 +=  dotmat(j, j-q, h) *  dotmat(j-k, j-q, l) * auxK;
              }
            }

            double auxK_outer = EpaK((t[j] - t[i]) / tau);


            temp += auxK_outer *temp1;
          }

          // Trace is scalar-equivalent for p=1
          tk(i, l) += temp;
        }
      }
    }

    // Normalization and global max update
    double Qtemp = arma::max(arma::abs(arma::vectorise(tk) / std::sqrt(n * tau * (sn - trim))));
    if (Qtemp > Q) Q = Qtemp;
  }

  return Q;
}


 //' @details
 //' difference based multiplier
 // [[Rcpp::export]]
 arma::mat Qdist_uni(arma::cube res,arma::cube dotmat, arma::vec t,  int trim,  int n, double tau,  int N, int B, int sn,  int L = 5){

   arma::mat R(n, B);
   R.randn();
   int q;
   arma::mat V;
   V =  Qmultiplierj_uni(tau, 0, trim, t, res, dotmat, sn);
   int rowV = V.n_rows;
   int colV = V.n_cols;
   // std::cout << rowV-2*L+1 << endl;
   int rowS = rowV - 2*L+1;
   arma::mat S(rowS, colV);
   arma::mat output(1,B);
   arma::mat boot(n - rowV - sn + 1, B),  tmp(N, B);
   int col1, col2;
   double frac;

   S.zeros();
   boot.zeros();
   output.zeros();

   for(int i = sn; i < L+sn; i++){
     S.row(sn) += V.row(i+L) - V.row(i); //L:2*L-1 - 0:L-1
   }
   for(int j = 1 + sn; j < rowS; j++){
     S.row(j) =   S.row(j - 1) + V.row(j + 2 * L - 1) + V.row(j-1) - 2*V.row(j + L - 1);
   }
   S/= sqrt(2*L);
   // std::cout<<max(vectorise(S))<<endl;
   for(int l = 0; l <= n - rowV - sn; l++){
     tmp.zeros();
     for(int j = sn; j < rowS; j++){
       col1 = l * N;
       col2 = l * N + N - 1;
       tmp += S.submat(j, col1, j, col2).t() * R.row(l+j);
     }
     tmp = abs(tmp);
     boot.row(l) = max(tmp); // column max
   }
   //  std::cout<<max(vectorise(boot))<<endl;
   frac = sn* (rowS) /2;
   //  std::cout<<frac<<endl;
   output = max(boot)/sqrt(frac);

   for(q = 1; q < N; q++){
     V =  Qmultiplierj_uni(tau, q, trim, t, res, dotmat, sn);
     S.zeros();
     boot.zeros();

     for(int i = sn; i < L+sn; i++){
       S.row(sn) += V.row(i+L) - V.row(i); //L:2*L-1 - 0:L-1
     }
     for(int j = 1 + sn; j < rowS; j++){
       S.row(j) =   S.row(j - 1) + V.row(j + 2 * L - 1) + V.row(j-1) - 2*V.row(j + L - 1);
     }
     S/= sqrt(2*L);
     // std::cout<<max(vectorise(S))<<endl;
     for(int l = 0; l <= n - rowV - sn; l++){
       tmp.zeros();
       for(int j = sn; j < rowS; j++){
         col1 = l * N;
         col2 = l * N + N - 1;
         tmp += S.submat(j, col1, j, col2).t() * R.row(l+j);
       }
       tmp = abs(tmp);
       boot.row(l) = max(tmp); // column max
     }
     //  std::cout<<max(vectorise(boot))<<endl;
     frac = sn* (rowS) /2;
     //  std::cout<<frac<<endl;
     output = max(output, max(boot)/sqrt(frac));
   }
   return(output);
 }

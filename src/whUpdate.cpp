/*
The core update procedure for Nonnegative Matrix Factorization via coordinate descent.
Author: Lanfeng Pan
Date: 2012-May
Using: RcppEigen

Compile with R CMD SHLIB whUpdate.cpp, then dyn.load("whupdate.so") in R.
or Build package with R CMD build bignmf then library(bignmf) in R
or 
  # g++ -I/usr/share/R/include -DNDEBUG -I/home/pkg/RcppEigen/include  -I/home/pkg/Rcpp/include   -fpic  -O3 -pipe  -g  -c myfile.cpp -o myfile.o
  # g++ -shared -o myfile.so myfile.o -L/home/pkg/Rcpp/lib -lRcpp -Wl,-rpath,/home/pkg/Rcpp/lib -L/usr/lib64/R/lib -lR
Then dyn.load in R
*/
  
  
  
#include <RcppEigen.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <time.h>

using namespace Eigen;
using namespace Rcpp;
typedef Map<MatrixXd> MapMatd;
//typedef Map<MatrixXi> MapMati;
//typedef Map<VectorXd> MapVecd;
/*Fast Matrix crossprod A' * A */
inline MatrixXd AtA_map(const MapMatd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
  .rankUpdate(A.adjoint());
}

/*Fast Matrix tcrossprod A* A' */
inline MatrixXd AAt_map(const MapMatd& A) {
  int n(A.rows());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
  .rankUpdate(A);
}

/* A'*A when A is MatrixXd type */
inline MatrixXd AtA(const MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
  .rankUpdate(A.adjoint());
}

/*A*A' when A is MatrixXd*/
inline MatrixXd AAt(const MatrixXd& A) {
  int n(A.rows());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
  .rankUpdate(A);
}


/* When updateing W and H, there is no lambda involved as glm via coordiante descent. Because I only need the solution when lambda=0. My experiments show it's even faster if I omit the steps from seq(max.lambda, 0, length=100)*/
/*Give H and V, update W.*/
void forW(const Eigen::MatrixXd& gram, Eigen::MatrixXd& gradw, Eigen::MatrixXd& w)
{
  Eigen::VectorXd delta(w.cols());
  Eigen::VectorXi sets(w.cols());
  
/*I jsut want sets = 1:w.cols(). Do I have to use loop in C++?*/
  for(int iset=0;iset<w.cols();iset++)
  {
    sets[iset]=iset;
  }
  int changetmp, changeidx;


  for(int irow =0; irow<w.rows(); irow++)
  {
    int iterations=0;
    while(iterations<2000)
    { 
      iterations++;
      delta.setZero(w.cols());
  /*Oops, here I just want rearange sets in random order, i.e. sets = sample.int(w.cols()). Updating in random order will speed up the convergence by several times */
	  srand((int)time(0));
      for(int k=0;k<w.cols();k++)
      {
        changeidx=rand() % w.cols();
        changetmp=sets[k];
        sets[k]=sets[changeidx];
        sets[changeidx]=changetmp;
      }
      
      for(int i =0;i<w.cols(); i++)
      {
        int icol=sets[i];

		/*If w(irow, icol) = w(irow, icol) - w(irow, icol)=0   */
		delta[icol]=std::min(w(irow, icol), gradw(irow, icol)/gram(icol, icol));
        if(std::abs(delta[icol]) > 1e-8)
        {
          gradw.row(irow) -= gram.row(icol) * delta[icol];
          w(irow, icol) -= delta[icol];         
        }
      }
      
      if(delta.squaredNorm() < 1e-8){
        break;
      }
    }
  }
  
}

void forH(const Eigen::MatrixXd& gram, Eigen::MatrixXd& gradh, Eigen::MatrixXd& h)
{
  Eigen::VectorXd delta(h.rows());
  Eigen::VectorXi sets(h.rows());
  for(int iset=0;iset<h.rows();iset++)
  {
    sets[iset]=iset;
  }
  int changetmp, changeidx;
  for(int icol =0; icol<h.cols(); icol++)
  {
    int iterations=0;
    while(iterations<2000)
    { 
      iterations++;
      delta.setZero(h.rows());
      srand((int)time(0));
      for(int k=0;k<h.rows();k++)
      {
        changeidx=rand() % h.rows();
        changetmp=sets[k];
        sets[k]=sets[changeidx];
        sets[changeidx]=changetmp;
      }
      
      for(int i =0;i<h.rows(); i++)
      {
        int irow=sets[i];
        delta[irow]=std::min(h(irow, icol), gradh(irow, icol)/gram(irow, irow));
        if(std::abs(delta[irow]) > 1e-8)
        {
          gradh.col(icol) -= gram.col(irow) * delta[irow];
          h(irow, icol) -= delta[irow];         
        }
      }
      if(delta.squaredNorm() < 1e-8){
        break;
      }
    }
  }
  
}

int wh(const MapMatd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& h, int& max_iter, double& tol)
{
  int n=v.rows();
  int m=v.cols();
  int r=w.cols();
  MatrixXd gram(r, r);
  MatrixXd gradh(r, m);
  MatrixXd gradw(n, r);
  
  double PgradeNorm=0;
  double initGrad=0;
  int iter;
  
  for(iter=0; iter<max_iter;iter++)
  {
    gram=AtA(w);
    gradh = gram * h;
    gradh.noalias() -= w.transpose() * v ;

/*Judge if the stop condition is satisfied*/
	if(iter==1)
    {
      initGrad= gradh.squaredNorm();
    }
    else if(iter > 1)
    {
      PgradeNorm=0;
      for(int icol =0; icol<gradh.cols();icol++)
      {
        for(int irow=0;irow<gradh.rows();irow++)
        {
          if(gradh(irow, icol)<0 || h(irow, icol)>0)
          {
            PgradeNorm += std::pow(gradh(irow, icol), 2);
          }
        }
      }
      if(PgradeNorm < initGrad * tol)
      {
        return iter;
      }
    }



	
    forH(gram, gradh, h);
    
    gram=AAt(h);
    gradw= w * gram;
    gradw.noalias() -= v * h.transpose();  
    forW(gram, gradw, w);
    
  }
  return iter;
}

/*This function maps data in R into Cpp*/
RcppExport SEXP whupdate(SEXP v_mat, SEXP w_mat, SEXP h_mat, SEXP max_iter_scalar, SEXP tol_scalar)
{
  BEGIN_RCPP
  
  const MapMatd v(as<MapMatd>(v_mat));
  MapMatd w_map(as<MapMatd>(w_mat));
  MapMatd h_map(as<MapMatd>(h_mat));
  int max_iter=Rcpp::as<int>(max_iter_scalar);
  double tol=Rcpp::as<double>(tol_scalar);
  tol=std::pow(tol, 2);
  MatrixXd w=w_map;
  MatrixXd h=h_map;
  int niter=0;
  niter=wh(v, w, h, max_iter, tol);
  
  return List::create(Named("W") = w, Named("H") = h, Named("iterations") = niter);
  
  END_RCPP
}


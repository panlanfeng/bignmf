/*
  # g++ -I/usr/share/R/include -DNDEBUG -I/home/pkg/RcppEigen/include  -I/home/pkg/Rcpp/include   -fpic  -O3 -pipe  -g  -c myfile.cpp -o myfile.o
  # g++ -shared -o myfile.so myfile.o -L/home/pkg/Rcpp/lib -lRcpp -Wl,-rpath,/home/pkg/Rcpp/lib -L/usr/lib64/R/lib -lR
  */
  
  
  
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;
inline MatrixXd AtA(const MapMatd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
  .rankUpdate(A.adjoint());
}

inline MatrixXd AAt(const MapMatd& A) {
  int n(A.rows());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
  .rankUpdate(A);
}




/*update beta and rx, one complete or incomptete circle*/
  int cdinnerloop(const Eigen::MatrixXd& gram, 
                  Eigen::VectorXd& rx, 
                  Eigen::VectorXd& beta, 
                  Eigen::VectorXi& sets
                  )
{
    int vk, mj,m;
    m=gram.cols();
    double threshold;
    //double sumofchange=0;
    int converge=0;
    Eigen::VectorXd change(m);
    
    /*Ramdomly disturbing the order of update  */
      int changetmp, changeidx;
    srand((int)time(0));
    for(int k=0;k<m;k++)
    {
      changeidx=rand() % m;
      changetmp=sets[k];
      sets[k]=sets[changeidx];
      sets[changeidx]=changetmp;
    }
    
    
    for(vk = 0; vk < m; vk++)
    {    
      mj=sets[vk];
      if(mj < 0 ) 
      {
        continue;
      }
      threshold = rx[mj] + beta[mj];
      
      if(threshold > 0)
      {
        change[mj] = beta[mj] - threshold;
        if(change[mj]!=0)
        {
            rx += gram.col(mj) * change[mj];
          /*beta change from zero to nonzero */
            /*update rx
            for( int newi= 0; newi < rx.size(); newi++)
            {
              rx[newi] += gram(newi, mj) * change[mj];
            }*/
          
          beta[mj]=threshold;  
        } 
      }
      else
      {
        change[mj]=beta[mj];
        if(change[mj]!=0)
        {
             rx += gram.col(mj) * change[mj];
          /*beta change from nonzero to zero 
            for( int newi= 0; newi < rx.size(); newi++)
            {
              rx[newi] += gram(newi, mj) * change[mj];
            }*/
          beta[mj]=0;
        }
      }
      
    }
    
    /*if the total change of beta is small, we think beta has been converged*/
      //sumofchange=change.squaredNorm()
     /* for(vk = 0; vk< m; vk++)
      {
        sumofchange = sumofchange + pow(change[vk], 2);
      }*/
    
    if(change.squaredNorm() < 1e-4)
    {
      converge=1;
    }
    return converge;
    
  }


/*one complete regression  */
Eigen::VectorXd cdouterloop(const Eigen::MatrixXd& gram, Eigen::VectorXd& rx, Eigen::VectorXd& beta)
{
    int converge=0;
    int m;
    m=gram.cols();
    //rx.noalias() -= gram * beta;
    Eigen::VectorXi allsets(m);
    for(int j=0;j<m;j++)
    {
      allsets[j]=j;
    }
    
    int iterations=0;
    /*iterate among all betas until converge */
      while(iterations<20000)
      { 
        iterations++;
        converge=cdinnerloop(gram, rx, beta, allsets);
        if(converge==1){
          return beta;
        }
      }
return beta;
}

RcppExport SEXP hupdate(SEXP v_mat, SEXP w_mat, SEXP h_mat)
{
  BEGIN_RCPP
  
  const MapMatd v(as<MapMatd>(v_mat));
  MapMatd w(as<MapMatd>(w_mat));
  MapMatd h_map(as<MapMatd>(h_mat));

  int n=v.rows();
  int m=v.cols();
  int r=w.cols();
  /*normalize gram wtv*/
  MatrixXd h=h_map;
  MatrixXd gram(r, r);
  //VectorXd xnorm(m);
 // MatrixXd norm_mat(m,m);
  MatrixXd wtv(r,m);
  gram=AtA(w);
  //xnorm=gram.diagonal().array().pow(-0.5);
  //norm_mat=xnorm * xnorm.transpose();
  //gram=norm_mat.array() * gram.array();
  wtv= w.transpose() * v;
  wtv.noalias() -=  gram * h;
  
  VectorXd tmpBeta(h.rows());
  VectorXd rx(wtv.rows());
  for(int icol=0;icol<h.cols();icol++)
  {
    tmpBeta=h.col(icol);
    rx=wtv.col(icol);
    h.col(icol)=cdouterloop(gram, rx, tmpBeta);
  }
  
  return wrap(h);
  
  END_RCPP
}

RcppExport SEXP wupdate(SEXP v_mat, SEXP w_mat, SEXP h_mat)
{
BEGIN_RCPP
  
  const MapMatd v(as<MapMatd>(v_mat));
  MapMatd h(as<MapMatd>(h_mat));
  MapMatd w_map(as<MapMatd>(w_mat));

  int n=v.rows();
  int m=v.cols();
  int r=h.rows();
  /*normalize gram wtv*/
  MatrixXd w=w_map;
  MatrixXd gram(r, r);
  //VectorXd xnorm(m);
 // MatrixXd norm_mat(m,m);
  MatrixXd vht(r,m);
  gram=AAt(h);
  //xnorm=gram.diagonal().array().pow(-0.5);
  //norm_mat=xnorm * xnorm.transpose();
  //gram=norm_mat.array() * gram.array();
  vht= v * h.transpose() ;
  vht.noalias() -= w * gram;  

  VectorXd tmpBeta(r);
  VectorXd rx(r); 
  for(int irow=0;irow<n;irow++)
  {
    tmpBeta=w.row(irow);
    rx=vht.row(irow);
    w.row(irow)=cdouterloop(gram, rx, tmpBeta);
  }
  
  return wrap(w);
END_RCPP
}




/* sparse matrix support*/
RcppExport SEXP sphupdate(SEXP v_mat, SEXP w_mat, SEXP h_mat)
{
  BEGIN_RCPP
  const MappedSparseMatrix<double> v(as<MappedSparseMatrix<double> >(v_mat));
  //const MapMatd v(as<MapMatd>(v_mat));
  MapMatd w(as<MapMatd>(w_mat));
  MapMatd h_map(as<MapMatd>(h_mat));

  int n=v.rows();
  int m=v.cols();
  int r=w.cols();
  /*normalize gram wtv*/
  MatrixXd h=h_map;
  MatrixXd gram(r, r);
  VectorXd xnorm(m);
  MatrixXd norm_mat(m,m);
  MatrixXd wtv(r,m);
  gram=AtA(w);
  xnorm=gram.diagonal().array().pow(-0.5);
  norm_mat=xnorm * xnorm.transpose();
  gram=norm_mat.array() * gram.array();
  wtv=xnorm.asDiagonal() * w.transpose() * v;
  
  VectorXd tmpBeta(h.rows());
  VectorXd rx(wtv.rows());
  for(int icol=0;icol<h.cols();icol++)
  {
    tmpBeta=h.col(icol);
    rx=wtv.col(icol);
    h.col(icol)=cdouterloop(gram, rx, tmpBeta);
  }
  
  return wrap(xnorm.asDiagonal() * h);
  
  END_RCPP
}

/*update sparse w matrix*/
RcppExport SEXP spwupdate(SEXP v_mat, SEXP w_mat, SEXP h_mat)
{
BEGIN_RCPP
  const MappedSparseMatrix<double> v(as<MappedSparseMatrix<double> >(v_mat));
  //const MapMatd v(as<MapMatd>(v_mat));
  MapMatd h(as<MapMatd>(h_mat));
  MapMatd w_map(as<MapMatd>(w_mat));

  int n=v.rows();
  int m=v.cols();
  int r=h.rows();
  /*normalize gram wtv*/
  MatrixXd w=w_map;
  MatrixXd gram(r, r);
  VectorXd xnorm(m);
  MatrixXd norm_mat(m,m);
  MatrixXd vht(r,m);
  gram=AAt(h);
  xnorm=gram.diagonal().array().pow(-0.5);
  norm_mat=xnorm * xnorm.transpose();
  gram=norm_mat.array() * gram.array();
  vht= v * h.transpose() * xnorm.asDiagonal();
  
  VectorXd tmpBeta(r);
  VectorXd rx(r); 
  for(int irow=0;irow<n;irow++)
  {
    tmpBeta=w.row(irow);
    rx=vht.row(irow);
    w.row(irow)=cdouterloop(gram, rx, tmpBeta);
  }
  
  return wrap(w * xnorm.asDiagonal());
END_RCPP
}


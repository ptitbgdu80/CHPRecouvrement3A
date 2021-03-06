#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>
#include "Sparse"
#include "Dense"
#include <vector>
#include <iostream>
#include <cmath>
#include "DataFile.h"
#include <string>
#include <fstream>

class Probleme
{
private:
  Eigen::SparseLU <Eigen::SparseMatrix<double> > _solver;
  Eigen::SparseMatrix<double> _Ap;
  Eigen::VectorXd _Up, _Bp, _CondBas, _CondHaut, _Utemps;
  int _NbLignes, _NbCol, _Np, _Me, _rec, _i1, _iN, _i1SansRec, _iNSansRec;
  int _choix;
  enum {stationnaire1, stationnaire2, instationnaire};
  int _formatSortie;
  enum {Paraview, Gnuplot, ParaviewEtGnuplot};
  double _Lx, _Ly, _Dx, _Dy, _Dt, _D, _tmax, _alpha, _beta, _C1, _C2, _C3, _t,_Epsilon;
  MPI_Status _Status;
  std::string _saveFolder;
public:
  Probleme(DataFile file);

  ~Probleme();

  void charge();

  void initializeSolver();

  double f(double x, double y, double t);

  double g(double x, double y);

  double h(double x, double y);

  void calculB();

  void communication();

  void SaveIteration();

  void TimeIteration();

  void PostProcessing();

  void CreationVtk();
};

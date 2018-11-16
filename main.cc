#include "fannyestlaplusjolie.h"
#include "DataFile.h"

int main(int argc, char * argv[])
{
  MPI_Status status;
  MPI_Init(&argc,&argv);
  int NbLignes, NbCol, rec;
  double Lx, Ly, Dt, tmax;

  NbLignes = 3;
  NbCol = 4;
  rec = 4;
  Lx = 1.;
  Ly = 1.;
  Dt = 0.1;
  tmax = 5.;
  

  // Probleme Pb1(DataFile);
  // Pb1.charge();
  // Pb1.initializeMatrix();

  MPI_Finalize();
  return 0;
}

#include "fannyestlaplusjolie.h"
#include "DataFile.h"

int main(int argc, char * argv[])
{
  MPI_Status status;
  MPI_Init(&argc,&argv);
  int NbLignes, NbCol, rec;
  double Lx, Ly, Dt, tmax;

  DataFile file(argv[1]);
  file.ReadDataFile();

  Probleme Pb1(file);
  Pb1.charge();
  Pb1.initializeMatrix();

  MPI_Finalize();
  return 0;
}

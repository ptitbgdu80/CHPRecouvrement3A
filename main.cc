#include "schwarz.h"
#include "DataFile.h"

int main(int argc, char * argv[])
{
  int Me;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Me);

  DataFile file(argv[1]);
  file.ReadDataFile();

  Probleme Pb1(file);
  Pb1.charge();
  Pb1.initializeSolver();

  Pb1.TimeIteration();

  Pb1.PostProcessing();

  MPI_Finalize();

  return 0;
}

#include "fannyestlaplusjolie.h"
#include "DataFile.h"

int main(int argc, char * argv[])
{
  MPI_Init(&argc,&argv);

  DataFile file(argv[1]);
  file.ReadDataFile();

  Probleme Pb1(file);
  Pb1.charge();
  Pb1.initializeSolver();

  MPI_Finalize();
  return 0;
}

#include "fannyestlaplusjolie.h"
#include "DataFile.h"

int main(int argc, char * argv[])
{
  int Me;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Me);
  double tstart, tend, time, t_total;

  DataFile file(argv[1]);
  file.ReadDataFile();

  tstart=MPI_Wtime();

  Probleme Pb1(file);
  Pb1.charge();
  Pb1.initializeSolver();

  Pb1.TimeIteration();

  tend=MPI_Wtime();
  t_total=tend-tstart;
  MPI_Allreduce(&t_total,&time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(Me==0)
  {
    printf("Total Elapsed Time = %f \n", time);
  }
  Pb1.PostProcessing();

  MPI_Finalize();

  return 0;
}

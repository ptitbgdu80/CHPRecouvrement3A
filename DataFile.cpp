#ifndef _DATA_FILE_CPP

#include "DataFile.h"

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name), _if_tmax(false), _if_Dt(false), _if_D(false),
_if_Lx(false), _if_Ly(false), _if_rec(false), _if_NbCol(false), _if_NbLignes(false), _if_CL(false)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &_Me);
}

void DataFile::ReadDataFile()
{
  ifstream data_file;
  data_file.open(_file_name);
  if (_Me == 0)
  {
    if (!data_file.is_open())
    {
      cout << "Unable to open file " << _file_name << endl;
      abort();
    }
    else
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Reading data file " << _file_name << endl;
    }
  }

  string file_line;

  while (!data_file.eof())
  {
    getline(data_file, file_line);

    if (file_line.find("Conditions aux limites") != std::string::npos)
    {
      data_file >> _alpha >> _beta; _if_CL = true;
    }

    if (file_line.find("choix") != std::string::npos)
    {
      data_file >> _choix; _if_tmax = true;
    }

    if (file_line.find("Dt") != std::string::npos)
    {
      data_file >> _Dt; _if_Dt = true;
    }

    if (file_line.find("tmax") != std::string::npos)
    {
      data_file >> _tmax; _if_tmax = true;
    }

    if (file_line.find("Nombre lignes") != std::string::npos)
    {
      data_file >> _NbLignes; _if_NbLignes = true;
    }

    if (file_line.find("Nombre colonnes") != std::string::npos)
    {
      data_file >> _NbCol; _if_NbCol = true;
    }

    if (file_line.find("recouvrement") != std::string::npos)
    {
      data_file >> _rec; _if_rec = true;
    }

    if (file_line.find("Dimension du domaine selon x") != std::string::npos)
    {
      data_file >> _Lx; _if_Lx = true;
    }

    if (file_line.find("Dimension du domaine selon y") != std::string::npos)
    {
      data_file >> _Ly; _if_Ly = true;
    }

    if (file_line.find("Coefficient de diffusion D") != std::string::npos)
    {
      data_file >> _D; _if_D = true;
    }
  }

  if (!_if_tmax)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (10) is used for tmax." << endl;
    }
    _tmax = 10.;
  }
  if (!_if_Dt)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (0.001) is used for Dt." << endl;
    }
    _Dt = 0.01;
  }
  if (!_if_choix)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default choice (stationnaire1) is used." << endl;
    }
    _choix = "stationnaire1";
  }
  if (!_if_CL)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default boundary condition is used." << endl;
    }
    _alpha = 0.;
    _beta = 1.;
  }
  if (!_if_rec)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value of the recovery is used for the velocity." << endl;
    }
    _rec = 3;
  }
  if (!_if_Lx)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value of the recovery is used for the velocity." << endl;
    }
    _Lx = 1.;
  }
  if (!_if_Ly)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value of the recovery is used for the velocity." << endl;
    }
    _Ly = 1.;
  }
  if (!_if_NbCol)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value of the recovery is used for the velocity." << endl;
    }
    _NbCol = 5;
  }
  if (!_if_NbLignes)
  {
    if (_Me == 0)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value of the recovery is used for the velocity." << endl;
    }
    _NbLignes= 5;
  }
}

#define _DATA_FILE_CPP
#endif

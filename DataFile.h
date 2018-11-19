#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <mpi.h>

// Définition de la classe

class DataFile {
private:
  std::string _file_name, _saveFolder;

  int _NbLignes, _NbCol, _rec, _Me;
  int _choix;
  enum{stationnaire1, stationnaire2, instationnaire};
  int _formatSortie;
  enum{Paraview, Gnuplot, ParaviewEtGnuplot};

  double _Lx, _Ly, _Dt, _tmax, _alpha, _beta, _D;

  bool _if_CL;
  bool _if_choix;
  bool _if_tmax;
  bool _if_Dt;
  bool _if_Lx;
  bool _if_Ly;
  bool _if_rec;
  bool _if_NbCol;
  bool _if_NbLignes;
  bool _if_D;
  bool _if_saveFolder;
  bool _if_formatSortie;

public: // Méthodes et opérateurs de la classe
  DataFile(std::string file_name);
  void ReadDataFile();
  std::string Get_file_name() const {return _file_name;};
  std::string Get_saveFolder() const {return _saveFolder;};
  double Get_tmax() const {return _tmax;};
  double Get_Dt() const {return _Dt;};
  double Get_Lx() const { return _Lx;};
  double Get_Ly() const { return _Ly;};
  double Get_D() const { return _D;};
  double Get_alpha() const { return _alpha;};
  double Get_beta() const { return _beta;};
  int Get_formatSortie() const {return _formatSortie;};
  int Get_choix() const {return _choix;};
  int Get_NbLignes() const { return _NbLignes;};
  int Get_NbCol() const { return _NbCol;};
  int Get_rec() const { return _rec;};

};

#define _DATA_FILE_H
#endif

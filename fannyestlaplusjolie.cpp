#include "fannyestlaplusjolie.h"
#include "DataFile.h"

Probleme::Probleme(DataFile file)
{
  _NbLignes = file.Get_NbLignes();
  _NbCol = file.Get_NbCol();
  _rec = file.Get_rec();
  _Lx = file.Get_Lx();
  _Ly = file.Get_Ly();
  _Dt = file.Get_Dt();
  _tmax = file.Get_tmax();
  _alpha = file.Get_alpha();
  _beta = file.Get_beta();
  _D = file.Get_D();
  MPI_Comm_rank(MPI_COMM_WORLD, &_Me);
  MPI_Comm_size(MPI_COMM_WORLD, &_Np);
  _Dx = _Lx/(_NbCol+1);
  _Dy = _Ly/(_NbLignes+1);
  _C1 = 1./_Dt + 2.*_D*(1./pow(_Dx,2)+1./pow(_Dy,2));
  _C2 = -_D/pow(_Dx,2);
  _C3 = -_D/pow(_Dy,2);
  _CondBas.Zero(_NbCol);
  _CondHaut.Zero(_NbCol);


  //Cette exception pose problème mais ça serait vraiment pas de bol...
  //Quand on essaye d'exprimer u_(i-1) en fonction de u_i avec la CL ça fait disparaître u_(i-1)
  if (_alpha-_beta*_Dy == 0)
  {
    std::cout << "***********************************************" << std::endl;
    std::cout << "***************alpha-beta*Dy = 0!**************" << std::endl;
    std::cout << "*******LA CONDTION DE BORD NE PERMET PAS*******" << std::endl;
    std::cout << "************DE RESOUDRE LE SYSTEME*************" << std::endl;
    std::cout << "***********************************************" << std::endl;
    std::cout << "alpha = " << _alpha << ", beta = " << _beta << ", Dy = " << _Dy << std::endl;
    exit(1);
  }
}

Probleme::~Probleme()
{}

void Probleme::charge()
{
  int q = (_NbLignes+2)/_Np;

  //rec = 2 recouvrement minimum (au moins deux lignes sont communes aux procs voisins)
  if (_rec < 2)
  {
    if (_Me == 0)
    {
      std::cout << "***********************************************" << std::endl;
      std::cout << "************RECOUVREMENT TROP PETIT************" << std::endl;
      std::cout << "***********************************************" << std::endl;
    }
    exit(1);
  }

  //rec = q recouvrement maximum (Au delà certaines lignes seront connues par 3 procs)
  if (_rec > q)
  {
    if (_Me == 0)
    {
      std::cout << "***********************************************" << std::endl;
      std::cout << "************RECOUVREMENT TROP GRAND************" << std::endl;
      std::cout << "***********************************************" << std::endl;
    }
    exit(1);
  }

  int r = _NbLignes+2 - q*_Np;
  int recHaut = _rec/2;
  int recBas = (_rec+1)/2;

  if (_Me < r)
  {
    _i1 = (q+1)*_Me;
    _iN = (q+1)*(_Me+1)-1;
  }
  else
  {
    _i1 = q*_Me + r;
    _iN = q*(_Me+1) + r - 1;
  }

  if (_Me != 0)
  {
    _i1 -= recBas;
  }
  if (_Me != _Np-1)
  {
    _iN += recHaut;
  }

  _Bp.Zero((_iN-_i1-1)*_NbCol);
  _Up.Zero((_iN-_i1-1)*_NbCol);
  _Ap.resize((_iN-_i1-1)*_NbCol,(_iN-_i1-1)*_NbCol);
}

void Probleme::initializeMatrix()
{
  std::vector<Eigen::Triplet<double>> liste_elem;

  for (int indice = 0; indice < (_iN-_i1-1)*_NbCol; indice++)
  {
    liste_elem.push_back({indice,indice,_C1});
  }

  for (int indice = 0; indice < (_iN-_i1-1)*_NbCol-1; indice++)
  {
    if ((indice+1)%_NbCol != 0)
    {
      liste_elem.push_back({indice,indice+1,_C2});
      liste_elem.push_back({indice+1,indice,_C2});
    }
  }

  for (int indice = 0 ; indice < (_iN-_i1-2)*_NbCol ; indice++)
  {
    liste_elem.push_back({indice,_NbCol+indice,_C3});
    liste_elem.push_back({_NbCol+indice,indice,_C3});
  }

  _Ap.setFromTriplets(liste_elem.begin(), liste_elem.end());

  //Bord bas
  for (int indice = 0; indice < _NbCol; indice ++)
  {
    //Le coeff vaut alors C1 + C3(alpha/(alpha-beta*Dy))
    _Ap.coeffRef(indice,indice) += _C3*_alpha/(_alpha-_beta*_Dy);
  }
  //Bord haut
  for (int indice = _NbCol*(_iN-_i1-2); indice < _NbCol*(_iN-_i1-1); indice ++)
  {
    //Le coeff vaut alors C1 + C3(alpha/(alpha+beta*Dy))
    _Ap.coeffRef(indice,indice) += _C3*_alpha/(_alpha+_beta*_Dy);
  }
}

double Probleme::f(double x, double y, double t)
{
  switch(_choix)
  {
    case stationnaire1:
    return 2.*(x - pow(x,2) + y - pow(y,2));
    break;

    case stationnaire2:
    return sin(x) + cos(y);
    break;

    case instationnaire:
    return exp(-pow((x-_Lx/2.),2))*exp(-pow((y-_Ly/2.),2))*cos(acos(-1)*t/2.);
    break;

    default:
    std::cout << "Problème avec la variable _choix" << std::endl;
    abort();
  }
}

double Probleme::g(double x, double y)
{
  switch(_choix)
  {
    case stationnaire1:
    return 0.;
    break;

    case stationnaire2:
    return sin(x) + cos(y);
    break;

    case instationnaire:
    return 0.;
    break;

    default:
    std::cout << "Problème avec la variable _choix" << std::endl;
    abort();
  }
}

double Probleme::h(double x, double y)
{
  switch(_choix)
  {
    case stationnaire1:
    return 0.;
    break;

    case stationnaire2:
    return sin(x) + cos(y);
    break;

    case instationnaire:
    return 1.;
    break;

    default:
    std::cout << "Problème avec la variable _choix" << std::endl;
    abort();
  }
}

void Probleme::calculB()
{
  for (int nl = _i1+1; nl < _iN; nl++)
  {
    for (int nc = 1; nc < _NbCol+1; nc++)
    {
      //Cas général B = f + u^n/Dt
      _Bp[(nl-_i1-1)*_NbCol + nc-1] = _Up[(nl-_i1-1)*_NbCol + nc-1]/_Dt + f(nc*_Dx,nl*_Dy,_t);

      if (nl == 1)
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1]+= _D*g(nc*_Dx,0)/pow(_Dy,2);
      }
      else if (nl == _NbLignes)
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1]+= _D*g(nc*_Dx,_Ly)/pow(_Dy,2);
      }
      else if (nl == _i1+1)
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1]+= _C3*_CondBas[nc-1]*_Dy/(_alpha-_beta*_Dy);
      }
      else if (nl == _iN-1)
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1]-= _C3*_CondHaut[nc-1]*_Dy/(_alpha+_beta*_Dy);
      }


      if (nc == 1)
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1]+= _D*h(0,nl*_Dy)/pow(_Dx,2);
      }
      else if (nc == _NbCol)
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1]+= _D*h(_Lx,nl*_Dy)/pow(_Dx,2);
      }
    }
  }
}

void Probleme::communication()
{
  Eigen::VectorXd tempHaut, tempBas; //Dans ces vecteurs on stocke la valeur de la condition de bord alpha dU + beta U qu'on enverra en haut et en bas
  tempBas.resize(_NbCol);
  tempHaut.resize(_NbCol);
  for (int nc = 1; nc <= _NbCol; nc++)
  {
    // K = alpha/Dy*(u_(i+1) - u_i) + beta*u_(i+1)
    tempHaut[nc-1] = _alpha/_Dy*(_Up[(_iN-_i1-_rec)*_NbCol] - _Up[(_iN-_i1-_rec-1)*_NbCol]) + _beta*_Up[(_iN-_i1-_rec)*_NbCol];
    // K = alpha/Dy*(u_i - u_(i-1)) + beta*u_(i-1)
    tempBas[nc-1] = _alpha/_Dy*(_Up[(_rec-1)*_NbCol] - _Up[(_rec-2)*_NbCol]) + _beta*_Up[(_rec-2)*_NbCol];
  }

  //Les procs pairs envoient puis reçoivent
  if (_Me%2 == 0)
  {
    if (_Me!=0)
    {
      MPI_Send(&tempHaut[0],_NbCol,MPI_DOUBLE,_Me-1,100,MPI_COMM_WORLD);
    }
    if (_Me!=_Np-1)
    {
      MPI_Send(&tempBas[0],_NbCol,MPI_DOUBLE,_Me+1,100,MPI_COMM_WORLD);
    }

    if (_Me!=0)
    {
      MPI_Recv(&_CondBas,_NbCol,MPI_DOUBLE,_Me-1,100,MPI_COMM_WORLD, &_Status);
    }
    if (_Me!=_Np-1)
    {
      MPI_Recv(&_CondHaut,_NbCol,MPI_DOUBLE,_Me+1,100,MPI_COMM_WORLD, &_Status);
    }
  }

  //Les procs impairs reçoivent puis envoient
  else
  {
    if (_Me!=0)
    {
      MPI_Recv(&_CondBas,_NbCol,MPI_DOUBLE,_Me-1,100,MPI_COMM_WORLD, &_Status);
    }
    if (_Me!=_Np-1)
    {
      MPI_Recv(&_CondHaut,_NbCol,MPI_DOUBLE,_Me+1,100,MPI_COMM_WORLD, &_Status);
    }

    if (_Me!=0)
    {
      MPI_Send(&tempHaut,_NbCol,MPI_DOUBLE,_Me-1,100,MPI_COMM_WORLD);
    }
    if (_Me!=_Np-1)
    {
      MPI_Send(&tempBas,_NbCol,MPI_DOUBLE,_Me+1,100,MPI_COMM_WORLD);
    }
  }
}

void Probleme::Rename(double t)
{
  std::string tn;
  int a,b,c;
  a = _Me/100;
  b =(_Me - 100*a)/10;
  c = _Me - 100*a -10*b;
  tn = std::to_string(a) + std::to_string(b) + std::to_string(c);
//  tn= 'bla';
  _savefile = "sol" + tn + "t" + std::to_string(t) + ".dat";
}

void Probleme::Save()
{
  int recHaut = _rec/2;
  int recBas = (_rec+1)/2;
  system(("mkdir -p ./" + std::to_string(_t)).c_str());
  std::ofstream mon_flux; // Contruit un objet "ofstream"
  mon_flux.open(_savefile, std::ios::out); // Ouvre un fichier appelé name_file
  if (_Me==0)
  {
    for (int i=_i1; i<_iN-_i1-recHaut-2; i++)
    {
      for(int j=0; j<_NbCol; j++)
      {
        mon_flux << j*_Dx << " " << i*_Dy << " " << _Up[i*_NbCol+j] << std::endl;
      }
    }
  }
  else if (_Me== _Np-1)
  {
    for (int i=recBas-1; i<_iN; i++)
    {
      for(int j=0; j<_NbCol; j++)
      {
        mon_flux << j*_Dx << " " << (_i1+i+1)*_Dy << " " << _Up[i*_NbCol+j] << std::endl;
      }
    }
  }
  else
  {
    for (int i=recBas-1; i<_iN-_i1-recHaut-2; i++)
    {
      for(int j=0; j<_NbCol; j++)
      {
        mon_flux << j*_Dx << " " << (_i1+i+1)*_Dy << " " << _Up[i*_NbCol+j] << std::endl;
      }
    }
  }
  mon_flux.close();
}

// void Probleme::vtk()
// {
//   mon_flux.open(name_file, std::ios::out); // Ouvre un fichier appelé name_file
//   mon_flux<<"# vtk DataFile Version 3.0\n"
//   <<"cell\n"
//   <<"ASCII\n"
//   <<"DATASET STRUCTURED_POINTS\n"
//   <<"DIMENSIONS "<< _numCols <<" "<<_numLines<<" 1\n"
//   <<"ORIGIN 0 0 0\n"
//   <<"SPACING 1.0 1.0 1.0\n"
//   <<"POINT_DATA "<<_numCols*_numLines<<"\n"
//   <<"SCALARS cell float\n"
//   <<"LOOKUP_TABLE default\n";
// }

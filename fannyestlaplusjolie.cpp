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
  _Epsilon = file.Get_Epsilon();
  _saveFolder = file.Get_saveFolder();
  _choix = file.Get_choix();
  _formatSortie = file.Get_formatSortie();
  MPI_Comm_rank(MPI_COMM_WORLD, &_Me);
  MPI_Comm_size(MPI_COMM_WORLD, &_Np);
  _Dx = _Lx/(_NbCol+1);
  _Dy = _Ly/(_NbLignes+1);
  if (_choix==stationnaire1 || _choix==stationnaire2)
  {
      _C1 = 2.*_D*(1./pow(_Dx,2)+1./pow(_Dy,2));
  }
  else
  {
      _C1 = 1./_Dt + 2.*_D*(1./pow(_Dx,2)+1./pow(_Dy,2));
  }
  _C2 = -_D/pow(_Dx,2);
  _C3 = -_D/pow(_Dy,2);
  _CondBas.setZero(_NbCol);
  _CondHaut.setZero(_NbCol);

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
    _i1SansRec = (q+1)*_Me;
    _iNSansRec = (q+1)*(_Me+1)-1;
  }
  else
  {
    _i1SansRec = q*_Me + r;
    _iNSansRec = q*(_Me+1) + r - 1;
  }

  _i1 = _i1SansRec;
  _iN = _iNSansRec;

  if (_Me != 0)
  {
    _i1 -= recBas;
  }
  if (_Me != _Np-1)
  {
    _iN += recHaut;
  }

  _Bp.setZero((_iN-_i1-1)*_NbCol);
  _Up.setZero((_iN-_i1-1)*_NbCol);
}

void Probleme::initializeSolver()
{
  if (_Me == 0)
  {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Initialisation des solveurs (création des matrices)" << std::endl;
  }

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


  _Ap.resize((_iN-_i1-1)*_NbCol,(_iN-_i1-1)*_NbCol);
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

  _solver.compute(_Ap);

  if (_Me == 0)
  {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Initialisation terminée" << std::endl;
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
    exit(1);
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
      if(_choix==stationnaire1 || _choix==stationnaire2)
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1] = f(nc*_Dx,nl*_Dy,_t);
      }
      else
      {
        _Bp[(nl-_i1-1)*_NbCol + nc-1] = _Up[(nl-_i1-1)*_NbCol + nc-1]/_Dt + f(nc*_Dx,nl*_Dy,_t);
      }
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
         _Bp[(nl-_i1-1)*_NbCol + nc-1]+= -_C3*_CondHaut[nc-1]*_Dy/(_alpha+_beta*_Dy);
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
  //  if (_Me==0)
  //  {
  //    std::cout << _Bp << std::endl;
  //  }
 }

void Probleme::communication()
{
  Eigen::VectorXd tempHaut, tempBas; //Dans ces vecteurs on stocke la valeur de la condition de bord alpha dU + beta U qu'on enverra en haut et en bas
  tempBas.setZero(_NbCol);
  tempHaut.setZero(_NbCol);

  for (int nc = 0; nc < _NbCol; nc++)
  {
      // K = alpha/Dy*(u_(i+1) - u_i) + beta*u_(i+1)
    if (_Me!=_Np-1)
    {
      tempHaut[nc] = (_alpha/_Dy)*(_Up[(_iN-_i1-_rec)*_NbCol+nc] - _Up[(_iN-_i1-_rec-1)*_NbCol+nc]) + _beta*_Up[(_iN-_i1-_rec)*_NbCol+nc];
    }
      // K = alpha/Dy*(u_i - u_(i-1)) + beta*u_(i-1)
    if (_Me!=0)
    {
      tempBas[nc] = (_alpha/_Dy)*(_Up[(_rec-1)*_NbCol+nc] - _Up[(_rec-2)*_NbCol+nc]) + _beta*_Up[(_rec-2)*_NbCol+nc];
    }
  }


  //Les procs pairs envoient puis reçoivent
  if (_Me%2 == 0)
  {

    if (_Me!=_Np-1)
    {
      MPI_Send(&tempHaut[0],_NbCol,MPI_DOUBLE,_Me+1,10*_Me,MPI_COMM_WORLD);
    }
    if (_Me!=0)
    {
      MPI_Send(&tempBas[0],_NbCol,MPI_DOUBLE,_Me-1,100*_Me,MPI_COMM_WORLD);
      //std::cout << _Me << " comm1 " << std::endl;
    }

    if (_Me!=0)
    {
      MPI_Recv(&_CondBas[0],_NbCol,MPI_DOUBLE,_Me-1,10*(_Me-1),MPI_COMM_WORLD, &_Status);
      //std::cout << _Me << " comm2 " << std::endl;
    }
    if (_Me!=_Np-1)
    {
      MPI_Recv(&_CondHaut[0],_NbCol,MPI_DOUBLE,_Me+1,100*(_Me+1),MPI_COMM_WORLD, &_Status);
      //std::cout << _Me << " comm2 " << std::endl;
    }
  }

  //Les procs impairs reçoivent puis envoient
  else
  {
    if (_Me!=0)
    {
      MPI_Recv(&_CondBas[0],_NbCol,MPI_DOUBLE,_Me-1,10*(_Me-1),MPI_COMM_WORLD, &_Status);
    }
    if (_Me!=_Np-1)
    {
      MPI_Recv(&_CondHaut[0],_NbCol,MPI_DOUBLE,_Me+1,100*(_Me+1),MPI_COMM_WORLD, &_Status);
      //std::cout << _Me << " comm3 " << std::endl;
    }

    if (_Me!=_Np-1)
    {
      MPI_Send(&tempHaut[0],_NbCol,MPI_DOUBLE,_Me+1,10*_Me,MPI_COMM_WORLD);
      //std::cout << _Me << " comm4 " << std::endl;
    }
    if (_Me!=0)
    {
      MPI_Send(&tempBas[0],_NbCol,MPI_DOUBLE,_Me-1,100*_Me,MPI_COMM_WORLD);
      //std::cout << _Me << " comm4 " << std::endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void Probleme::TimeIteration()
{
  Eigen::VectorXd Up_Avant;
  double erreurloc,erreur;
  int iter=0;

  if (_Me == 0)
  {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Boucle temporelle de résolution du problème" << std::endl;
  }
  _t=0.;
  if (_Me == 0)
  {
    system(("mkdir -p " + _saveFolder).c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);
  while (_t<_tmax)
  {
    iter=0.;
    SaveIteration();
    _t += _Dt;
    erreur = _Epsilon+1.;
    while(erreur > _Epsilon)
    {
      iter=iter+1;
      communication();
      calculB();
      Up_Avant=_Up;
      _Up=_solver.solve(_Bp);
      erreurloc=(_Up-Up_Avant).norm();
      MPI_Allreduce(&erreurloc,&erreur,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      if(_Me==0)
      {
        std::cout << erreur <<std::endl;
      }
    }
    std::cout << "_t = " << _t << std::endl;
  }
  if (_Me == 0)
  {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Résolution effectuée avec succès" << std::endl;
  }
  SaveIteration();
  std::cout << "iter =" << iter << std::endl;
}

void Probleme::SaveIteration()
{
  int a,b,c;
  a = _Me/100;
  b =(_Me - 100*a)/10;
  c = _Me - 100*a -10*b;
  std::string procID = std::to_string(a) + std::to_string(b) + std::to_string(c);
  std::string savefile = "t" + std::to_string(_t) + "p" + procID + ".dat";

  std::ofstream mon_flux ; // Contruit un objet "ofstream"
  mon_flux.open( _saveFolder + "/" + savefile, std::ios::out); // Ouvre un fichier appelé name_file
  if (_Me==0 && _Me==_Np-1)
  {
    //On ne prend ni la ligne 0 ni la ligne iN
    for (int nl=1; nl < _iNSansRec; nl++)
    {
      for (int nc = 1; nc < _NbCol+1; nc++)
      {
        //Le premier élément de Up à prendre est celui d'indice 0
        mon_flux << nc*_Dx << " " << nl*_Dy << " " << _Up[(nl-1)*_NbCol+nc-1] << std::endl;
      }
    }
  }
  else if (_Me==0)
  {
    //On ne prend pas la ligne 0 (bord bas du domaine), on prend cependant la ligne iNSansRec
    // std::cout << "iNSansRec = " << _iNSansRec << std::endl;
    for (int nl=1; nl < _iNSansRec+1; nl++)
    {
      for(int nc = 1; nc < _NbCol+1; nc++)
      {
        // std::cout << nl << " " << nc << std::endl;
        //Le premier élément de Up à prendre est celui d'indice 0
        mon_flux << nc*_Dx << " " << nl*_Dy << " " << _Up[(nl-1)*_NbCol+nc-1] << std::endl;
      }
    }
  }
  else if (_Me== _Np-1)
  {
    //On ne prend pas la dernière ligne (iNSansRec = iN correspond au bord haut du domaine), on prend la ligne i1SansRec
    for (int nl = _i1SansRec; nl < _iNSansRec; nl++)
    {
      for(int nc = 1; nc < _NbCol+1; nc++)
      {
        //Le premier élément de Up à prendre est celui d'indice (i1SansRec - i1 - 1)*NbCol
        mon_flux << nc*_Dx << " " << nl*_Dy << " " << _Up[(nl-_i1-1)*_NbCol + nc-1] << std::endl;
      }
    }
  }
  else
  {

    //On prend les deux lignes i1SansRec et iNSansRec
    for (int nl = _i1SansRec; nl < _iNSansRec+1; nl++)
    {
      for(int nc = 1; nc < _NbCol+1; nc++)
      {
        //Le premier élément de Up à prendre est celui d'indice (i1SansRec - i1 - 1)*NbCol
        mon_flux << nc*_Dx << " " << nl*_Dy << " " << _Up[(nl-_i1-1)*_NbCol + nc-1] << std::endl;
      }
    }
  }
  mon_flux.close();
}

void Probleme::PostProcessing()
{
  if (_Me == 0)
  {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Début du postprocessing" << std::endl;
  }
  switch (_formatSortie)
  {
    case Paraview:
    CreationVtk();
    break;

    // case Gnuplot:
    //
    // break;

    // case ParaviewEtGnuplot:
    // CreationVtk();
    //
    // break;

    default:
    std::cout << "Problème avec le format de sortie" << std::endl;
    exit(1);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (_Me==0)
  {
    system(("rm -rf " + _saveFolder).c_str());
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Fin du postprocessing" << std::endl;
  }

}

void Probleme::CreationVtk()
{
  if (_Me == 0)
  {
    system(("rm -rf " + _saveFolder + "Paraview").c_str());
    system(("mkdir -p " + _saveFolder + "Paraview").c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //permet de démarrer avec un _t "multiple" de _Dt
  double tmin = floor(_Me*_tmax/(_Np*_Dt))*_Dt;
  double tmax = floor((_Me+1)*_tmax/(_Np*_Dt))*_Dt;
  std::cout << tmin << " " << tmax << std::endl;
  if (_Me == _Np-1)
  {
    tmax = _tmax + _Dt/2;
  }

  for (_t = tmin; _t < tmax; _t += _Dt)
  {
    //+0.1 pour éviter les erreurs d'arrondis, sachant que t/Dt est censé être entier
    int it = floor(_t/_Dt+0.1);

    std::ofstream mon_flux;
    mon_flux.open(_saveFolder + "Paraview/it" + std::to_string(it) + ".vtk", std::ios::out);
    mon_flux << "# vtk DataFile Version 3.0\n"
    <<"cell\n"
    <<"ASCII\n"
    <<"DATASET STRUCTURED_POINTS\n"
    <<"DIMENSIONS " << _NbCol << " " << _NbLignes << " 1\n"
    <<"ORIGIN 0 0 0\n"
    <<"SPACING  " << _Dx << " " << _Dy <<  "1.0\n"
    <<"POINT_DATA " << _NbCol*_NbLignes << "\n"
    <<"SCALARS cell float\n"
    <<"LOOKUP_TABLE default";

    for (int i = 0; i < _Np; i++)
    {
      int a,b,c;
      a = i/100;
      b =(i - 100*a)/10;
      c = i - 100*a -10*b;
      std::string procID = std::to_string(a) + std::to_string(b) + std::to_string(c);
      std::string savefile = "t" + std::to_string(_t) + "p" + procID + ".dat";
      std::ifstream fichier(_saveFolder + "/" + savefile);

      std::string x, y, newy, valeur;
      fichier >> x >> newy >> valeur;
      y = "";

      while (!fichier.eof())
      {
        if (newy == y)
        {
          mon_flux << valeur << " ";
        }
        else
        {
          mon_flux << std::endl << valeur << " ";
        }
        y = newy;
        fichier >> x >> newy >> valeur;
      }
    }
    mon_flux.close();
  }
}

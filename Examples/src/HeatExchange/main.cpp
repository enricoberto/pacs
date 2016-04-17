#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
#include "norm.hpp"
#include "matrid.hpp"

/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

using namespace std; // avoid std::

//! main program
int main(int argc, char** argv)
{
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  const int&    itermax= param.itermax;   //max number of iteration for Gauss-Siedel
  const double& toler=param.toler;   // Tolerance for stopping criterion
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto& M=param.M; // Number of grid elements
  const auto& name_out=param.n_o;
  const auto& mode=param.m; // Mode
  const auto& criterion=param.cr; // Criterion
  const auto& dt=param.dt; // Delta t
  const auto& iter_t=param.itt; // Iter tempo
  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1);
  std::vector<double> v(M+1);
  
  // Gauss Siedel is initialised with a linear variation
  // of T
  
  for(unsigned int m=0;m <= M;++m)
     theta[m]=(1.-m*h)*(To-Te)/Te;
 
  // Thomas
  // Inserisco la soluzione nella quarta colonna dei risultati.
	vector<double> diagi(M-1,-1);
	vector<double> diag(M, (2+h*h*act));
	vector<double> diags(M-1, -1);
	diag[M-1]=1;
	vector<double> b(M, 0);
	vector<double> u(M, 0);
	b[0]= theta[0];

	matrid A(diagi, diag, diags);
	A.thomas_algorithm(b, u);

	u.insert(u.begin(), theta[0]);


	
  // Gauss-Seidel
  // epsilon=||x^{k+1}-x^{k}||
  // Stopping criteria epsilon<=toler
  
	cout << "Il criterio di arresto e' "; 
  int iter=0;
  double xnew, epsilon;
     do
       { epsilon=0.;

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   v[m] = xnew - theta[m];
	   //epsilon += (xnew-theta[m])*(xnew-theta[m]);
	   theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1];
	 v[M] = xnew - theta[M];
	 //epsilon += (xnew-theta[M])*(xnew-theta[M]);
	 theta[M]=  xnew; 

	 iter=iter+1;    

 	 // compute norms
	switch(criterion){
		case 1:
			if(iter==1) cout << "Rn norm" << endl;
			for (unsigned int m = 1; m <=M; ++m) {
				epsilon += v[m]*v[m];
			}
			break;
		case 2:
			if(iter==1) cout << "L2 norm" << endl;
			epsilon = norm_L2(v, h);
			epsilon = norm_L2(v, h);
			break;
		case 3:
			if(iter==1) cout << "H1 norm" << endl;
			epsilon = norm_H1(v, h);
			break;
		default:
			if(iter==1) cout << "The criterion is wrong" << endl;
	 
	 }

       }while((sqrt(epsilon) > toler) && (iter < itermax) );

    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }

 // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));


	//TIME EVOLUTION
	
	// (I + dt*A) theta n+1 = b*dt + theta n
	// inizializzazione matrice A in cui aggiugo una prima linea con 1 in pos 0
	// per imporre la condizione al bordo
	vector<double> diagTi(M,-1);
	vector<double> diagT(M+1, (2+h*h*act));
	vector<double> diagTs(M, -1);
	diagT[M]=1;
	diagTi[0]=0;
	diagT[0]=1;
	diagTs[0]=0;

	matrid AT(diagTi, diagT, diagTs);


	// inizializzazione termine noto e vettore soluzioni
	vector<double> bT(M+1, 0);
	bT[0]= theta[0]; //b=b*dt
	bT[1]= theta[0]; //b=b*dt

	vector<double> theta_n(M+1, (To-Te)/Te); //cond iniziale: temp esterna
	vector<double> theta_n1(M+1, 0);

	AT.time_evolution(theta_n, theta_n1,bT, dt, iter_t);

	//PLOT

     // writing results with format
     // x_i u_h(x_i) u(x_i) uthomas(x_i) u_time(x_i) and lauch gnuplot 

   Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);
     std::vector<double> thomas(M+1);
     std::vector<double> time(M+1);
	if(mode==0 || mode==2){
		 const char *name_punt = name_out.c_str();

		cout<<"Result file: " << name_punt <<endl;
		ofstream f(name_punt);
    
	   	for(int m = 0; m<= M; m++)
		 {
	 // \t writes a tab 
		    f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<"\t"<<Te*(1+ u[m]) <<"\t"<<Te*(1+ theta_n1[m]) <<endl;
		 }
    f.close();
	}
     // Using temporary files (another nice use of tie)
    if(mode==1 || mode==2){
	   	for(int m = 0; m<= M; m++){
	 // An example of use of tie and tuples!
         
			std::tie(coor[m],sol[m],exact[m],thomas[m],time[m])=
			std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m],Te*(1.+ u[m]),Te*(1.+ theta_n1[m]));
		}
		gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex',"<< gp.file1d(std::tie(coor,thomas))<<
       "w l title 'uthom'," << gp.file1d(std::tie(coor,time))<<
       "w l title 'u after T seconds'" <<std::endl;
	}
     return status;
}

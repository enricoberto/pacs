#include "matrid.hpp"
//#include <cassert>
//#include <algorithm>
//#include <ctime>
//#include <Eigen/Dense>

/*
matrid matrid::transpose () const
{
	matrid mat(*this);
	
	std::vector<double> temp = this->get_a ();
	this->set_a(this.get_c);
	this->set_c(temp);

	return mat;
}

*/

/*
matrid operator* (const matrid A, const double k){
	std::cout << "prodotto" << std::endl;	
	std::vector<double> vec1(A.get_a());
	for (unsigned int i = 0; i < vec1.size(); ++i){vec1[i]*=k;}
	std::vector<double> vec2(A.get_b());
	for (unsigned int i = 0; i < vec2.size(); ++i){vec2[i]*=k;}
	std::vector<double> vec3(A.get_c());
	for (unsigned int i = 0; i < vec3.size(); ++i){vec3[i]*=k;}

	matrid B(vec1, vec2, vec3);
	 return B;
 }
*/
void matrid::scalar(const double k){
	for (unsigned int i = 0; i < N-1; ++i){
		this->mod_a()[i] *=k;
		this->mod_b()[i] *=k;
		this->mod_c()[i] *=k;
	}	
		this->mod_b()[N-1] *=k;
}
/*
 matrid& operator=(matrid &rhs){
	this->mod_a = rhs->get_a;
	this->mod_b = rhs->get_b;
	this->mod_c = rhs->get_c;

	return *this;
 }
 matrid operator=(matrid rhs){
	this.mod_a = rhs->get_a;
	this.mod_b = rhs->get_b;
	this->mod_c = rhs->get_c;

	return this;
 }
*/

 void matrid::set_gamma(){
     // creo vettori aT bT cT per uniformarmi alla notazione del
	 // Quarteroni-Sacco

	 std::vector<double>  bT = a;
	 std::vector<double>  aT = b;
	 std::vector<double>  cT = c;
	 
	 bT.insert(bT.begin(), 0);
	 cT.push_back(0);
	 
	 gamma[0]=1/aT[0];
	for (unsigned int i = 1; i < N; ++i) {
		gamma[i]=1/(aT[i]-bT[i]*gamma[i-1]*cT[i-1]);
	} 
	gamma_set=true;
 }

 void matrid::thomas_algorithm(const std::vector<double>& f,   std::vector<double>& x){
     // creo vettori aT bT cT per uniformarmi alla notazione del
	 // Quarteroni-Sacco
	 std::vector<double>  bT = a;
	 std::vector<double>  aT = b;
	 std::vector<double>  cT = c;
	 size_t n = f.size();
	
	 std::vector<double>  gammaT(n, 0);
	 std::vector<double>  y(n, 0);
    
	 bT.insert(bT.begin(), 0);
	 cT.push_back(0);
	
	//set gamma
	if(!gamma_set) this->set_gamma();
	y[0]=gamma[0]*f[0];
	
	for (unsigned int i = 1; i < n; ++i) {
		y[i]=gamma[i]*(f[i]-bT[i]*y[i-1]);		
	} 
	x[n-1]=y[n-1];
	for (unsigned int i = n-1; i-- >0;) {
		x[i]=y[i]-gamma[i]*cT[i]*x[i+1];
	}

 }


void matrid::time_evolution( std::vector<double>& theta_n,  std::vector<double>& theta_n1,  std::vector<double>& b, double dt, int iter_t){
	//size
	int n(theta_n.size() );

	// mat = mat*dt + I
	this->scalar(dt);
	for (unsigned int i = 1; i < n; ++i) this->mod_b()[i] +=1;
	// mat[0,0]=1
	this->mod_b()[0]=1;

	//termine noto
	std::vector<double> rhsT(n, 0);

	//Iterazioi nel tempo. Sistema lineare da risolvere: AT * theta_{n+1}=rhsT
	for (unsigned int it_t = 0; it_t < iter_t; ++it_t) {
		 //termine noto : b + theta;	
		for (unsigned int i = 0; i < n; ++i) rhsT[i] = b[i]*dt + theta_n[i];  
		rhsT[0]=b[0];

		// risoluzione sistema lineare
		this->thomas_algorithm(rhsT, theta_n1);
		theta_n = theta_n1;
	
		//cout
	//	for (unsigned int i = 0; i < n; ++i) std::cout << theta_n1[i] << " ";
	//	std::cout << std::endl;
	}

}





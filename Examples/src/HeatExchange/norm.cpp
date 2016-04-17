#include "norm.hpp"

double norm_L2(std::vector<double> v, double h){
	double sum(0.);
	std::vector<double> w(v.size(), 1.);
	w[0]=0.5; w[w.size()-1]=0.5;
	for (unsigned int i = 0; i < v.size(); ++i) {
		sum += w[i]*v[i]*v[i]*h;
	}
	return sum;
}


double norm_H1(std::vector<double> v, double h){

	double L2norm = norm_L2(v, h);

	std::vector<double> vp(v.size());
	
	vp[0]=(v[1]-v[0])/h;
	for (unsigned int i = 1; i < vp.size()-1; ++i) {
		vp[i]=(v[i+1]-v[i-1])/(2*h);	
	}
	vp[vp.size()-1]=(v[v.size()-1]-v[v.size()-2])/h;

	double L2normp = norm_L2(vp, h);

	return  (L2normp + L2norm);

}

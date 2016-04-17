#ifndef HAVE_MATRID_H
#define HAVE_MATRID_H

#include <vector>
#include <iostream>

class
matrid
{

private :
  
  const unsigned int N;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
  
/*  inline
  unsigned int
  sub2ind (const unsigned int ir,
           const unsigned int jc) const
  { return (ir + jc * rows); };
*/

public :

  matrid (unsigned int size)
    : N (size)
  {
	  a.resize (size-1, 0.0);
	  b.resize (size, 0.0);
	  c.resize (size-1, 0.0);
  }  	
  
  matrid (std::vector<double> &vec1, 
		  std::vector<double> &vec2,
		  std::vector<double> &vec3):  N(vec2.size()), a(vec1), b(vec2), c(vec3){}


  matrid (matrid const &) = default;  
  
  unsigned int
  get_size () const { return N; }

  double &
  operator() (unsigned int irow, unsigned int jcol)
  {
	  if(irow==(jcol+1)) return a[jcol];
	  else if(irow==jcol) return b[irow];
	  else if(irow==(jcol-1)) return c[irow];
	  else std::cerr << "NOOOOOOOOOOOO" << std::endl;

  };

  const double &
  operator()  (unsigned int irow, unsigned int jcol) const
  {
	  if(irow==(jcol+1)) return a[jcol];
	  else if(irow==jcol) return b[irow];
	  else if(irow==(jcol-1)) return c[irow];
  };

  
  double  *   mod_a () { return &(a[0]); };
  double  *   mod_b () { return &(b[0]); };
  double  *   mod_c () { return &(c[0]); };


  std::vector<double>   get_a () { return a; };
  const std::vector<double>  get_a () const { return a; };
  std::vector<double>  get_b () { return b; };
  const std::vector<double>    get_b () const { return b; };
  std::vector<double>    get_c () { return c; };
  const std::vector<double>    get_c () const { return c; };
  void  set_a (std::vector<double> vec) { a=vec;};
  void  set_b (std::vector<double> vec) { b=vec;};
  void  set_c (std::vector<double> vec) { c=vec;};
 
  
  void thomas_algorithm( const std::vector<double>& f,  std::vector<double>& x); 
  
void time_evolution( std::vector<double>& theta_n,  std::vector<double>& theta_n1, std::vector<double>& b, double dt, int iter_t);

// matrid& operator=(matrid &rhs);
 
//matrid& operator* (const double k);

void scalar(const double k);
};

#endif

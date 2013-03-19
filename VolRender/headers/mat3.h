// mat3.h: interface for the mat3 class.
//
// Defines a 3x3 matrix
// main operations are: multiplication, addition, multiplication with vector
//
// By JMD 9/8/06
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MAT3_H__1544EA11_662A_41FD_8BD3_D7311F73A131__INCLUDED_)
#define AFX_MAT3_H__1544EA11_662A_41FD_8BD3_D7311F73A131__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vec3.h"


template <class T> class mat3 // : public hviAlgebra<mat3<T>,T>  
{
protected:
	vec3<T> i,j,k;

public:
	// Constructor: defines an identity matrix
	mat3<T>():i(vec3<T>(1,0,0)),j(vec3<T>(0,1,0)),k(vec3<T>(0,0,1)) {   }
	// defines a matrix by three 3D vectors each representing a column 
	mat3<T>(const vec3<T> &a, const vec3<T> &b, const vec3<T> &c):i(a),j(b),k(c) { }

	// defines a matrix with a single vector u as : M=uu^t (vector u times its transposed u^t) 
	mat3<T>(const vec3<T> &u) 
	{
		vec3<T> v(u);
		i.scale(u, v.X());
		j.scale(u, v.Y());
		k.scale(u, v.Z());
	}

	// Selectors : vectors corresponding to columns
	vec3<T> I() const { return i; }
	vec3<T> J() const { return j; }
	vec3<T> K() const { return k; }

	// multiply matrix by vector, result is vector
	vec3<T> mult(const vec3<T> &v) const
	{
		return vec3<T> (	i.X()*v.X()+j.X()*v.Y()+k.X()*v.Z() ,
							i.Y()*v.X()+j.Y()*v.Y()+k.Y()*v.Z() ,
							i.Z()*v.X()+j.Z()*v.Y()+k.Z()*v.Z() );
	}
	vec3<T> operator*(const vec3<T> &v) const
	{
		return vec3<T> (	i.X()*v.X()+j.X()*v.Y()+k.X()*v.Z() ,
							i.Y()*v.X()+j.Y()*v.Y()+k.Y()*v.Z() ,
							i.Z()*v.X()+j.Z()*v.Y()+k.Z()*v.Z() );
	}

	// multiply by a scalar value all components
	void scale(T x)
	{
		i.scale(x);
		j.scale(x);
		k.scale(x);
	}
	void operator*=(T x)
	{
		i.scale(x);
		j.scale(x);
		k.scale(x);
	}
	void operator/=(T x)
	{
		i/=x;
		j/=x;
		k/=x;
	}
	void scale(const mat3<T> &m, T x)
	{
		i.scale(m.i, x);
		j.scale(m.j, x);
		k.scale(m.k, x);
	}

	// multiply two matrices
	void mult(const mat3<T> &a, const mat3<T> &b)
	{
		mat3<T> r;
		r.i=vec3<T>(	a.i.X()*b.i.X()+a.j.X()*b.i.Y()+a.k.X()*b.i.Z(),
						a.i.Y()*b.i.X()+a.j.Y()*b.i.Y()+a.k.Y()*b.i.Z(),
						a.i.Z()*b.i.X()+a.j.Z()*b.i.Y()+a.k.Z()*b.i.Z() );

		r.j=vec3<T>(	a.i.X()*b.j.X()+a.j.X()*b.j.Y()+a.k.X()*b.j.Z(),
						a.i.Y()*b.j.X()+a.j.Y()*b.j.Y()+a.k.Y()*b.j.Z(),
						a.i.Z()*b.j.X()+a.j.Z()*b.j.Y()+a.k.Z()*b.j.Z() );

		r.k=vec3<T>(	a.i.X()*b.k.X()+a.j.X()*b.k.Y()+a.k.X()*b.k.Z(),
						a.i.Y()*b.k.X()+a.j.Y()*b.k.Y()+a.k.Y()*b.k.Z(),
						a.i.Z()*b.k.X()+a.j.Z()*b.k.Y()+a.k.Z()*b.k.Z() );
		*this = r;
	}

	// multiply two matrices
	mat3<T> operator*(const mat3<T> &b) const
	{
		mat3<T> r;

		r.i=vec3<T>(	i.X()*b.i.X()+j.X()*b.i.Y()+k.X()*b.i.Z(),
						i.Y()*b.i.X()+j.Y()*b.i.Y()+k.Y()*b.i.Z(),
						i.Z()*b.i.X()+j.Z()*b.i.Y()+k.Z()*b.i.Z() );

		r.j=vec3<T>(	i.X()*b.j.X()+j.X()*b.j.Y()+k.X()*b.j.Z(),
						i.Y()*b.j.X()+j.Y()*b.j.Y()+k.Y()*b.j.Z(),
						i.Z()*b.j.X()+j.Z()*b.j.Y()+k.Z()*b.j.Z() );

		r.k=vec3<T>(	i.X()*b.k.X()+j.X()*b.k.Y()+k.X()*b.k.Z(),
						i.Y()*b.k.X()+j.Y()*b.k.Y()+k.Y()*b.k.Z(),
						i.Z()*b.k.X()+j.Z()*b.k.Y()+k.Z()*b.k.Z() );
		return r;
	}

	// multiply two matrices
	void operator*=(const mat3<T> &b)
	{
		mat3<T> r;

		r.i=vec3<T>(	i.X()*b.i.X()+j.X()*b.i.Y()+k.X()*b.i.Z(),
						i.Y()*b.i.X()+j.Y()*b.i.Y()+k.Y()*b.i.Z(),
						i.Z()*b.i.X()+j.Z()*b.i.Y()+k.Z()*b.i.Z() );

		r.j=vec3<T>(	i.X()*b.j.X()+j.X()*b.j.Y()+k.X()*b.j.Z(),
						i.Y()*b.j.X()+j.Y()*b.j.Y()+k.Y()*b.j.Z(),
						i.Z()*b.j.X()+j.Z()*b.j.Y()+k.Z()*b.j.Z() );

		r.k=vec3<T>(	i.X()*b.k.X()+j.X()*b.k.Y()+k.X()*b.k.Z(),
						i.Y()*b.k.X()+j.Y()*b.k.Y()+k.Y()*b.k.Z(),
						i.Z()*b.k.X()+j.Z()*b.k.Y()+k.Z()*b.k.Z() );
		*this= r;
	}

	// divide two matrices (multiply with inverse)
	void div(const mat3<T> &a, const mat3<T> &b)
	{
		mat3<T> r(b);
		T d = r.det();
		if (d==T(0)) { printf("cannot divide by matrice!"); exit(1); }
		r.inverse(r, d);
		mult(a, r);
	}
	void operator/=(const mat3<T> &b)
	{
		mat3<T> r(b);
		T d = r.det();
		if (d==T(0)) { printf("cannot divide by matrice!"); exit(1); }
		r.inverse(r, d);
		mult(*this, r);
	}
	mat3<T> operator/(const mat3<T> &b) const
	{
		mat3<T> r(b);
		T d = r.det();
		if (d==T(0)) { printf("cannot divide by matrice!"); exit(1);  }
		r.inverse(r, d);
		r.mult(*this, r);
		return r;
	}

	// add two matrices
	void add(const mat3<T> &a, const mat3<T> &b)
	{
		i.add(a.i, b.i);
		j.add(a.j, b.j);
		k.add(a.k, b.k);
	}
	mat3<T> operator+(const mat3<T> &b) const
	{
		mat3<T> r;

		r.i=i+b.i;
		r.j=j+b.j;
		r.k=k+b.k;
		return r;
	}
	void operator+=(const mat3<T> &b)
	{
		i+=b.i;
		j+=b.j;
		k+=b.k;
	}

	// sub two matrices
	void sub(const mat3<T> &a, const mat3<T> &b)
	{
		i.sub(a.i, b.i);
		j.sub(a.j, b.j);
		k.sub(a.k, b.k);
	}
	mat3<T> operator-(const mat3<T> &b) const
	{
		mat3<T> r;

		r.i=i-b.i;
		r.j=j-b.j;
		r.k=k-b.k;
		return r;
	}
	void operator-=(const mat3<T> &b)
	{
		i-=b.i;
		j-=b.j;
		k-=b.k;
	}


	// compute determinant (a scalar value)
	T det() const
	{
		return	  i.X()*(j.Y()*k.Z()-k.Y()*j.Z()) 
				- i.Y()*(j.X()*k.Z()-k.X()*j.Z())
				+ i.Z()*(j.X()*k.Y()-k.X()*j.Y());
	}

	// Inverse the matrix a, 
	// works only if det/=0, det must be the determinant of a
	void inverse(const mat3<T> &a, T det)
	{
		if (det==T(0)) { printf("cannot inverse with nul det!"); exit(1);  }
		vec3<T> ai(a.i), aj(a.j), ak(a.k);
		i = vec3<T>(	  aj.Y()*ak.Z()-ak.Y()*aj.Z(),
						-(aj.X()*ak.Z()-ak.X()*aj.Z()),
						  aj.X()*ak.Y()-ak.X()*aj.Y()  );

		j = vec3<T>(	-(ai.Y()*ak.Z()-ak.Y()*ai.Z()),
						  ai.X()*ak.Z()-ak.X()*ai.Z(),
						-(ai.X()*ak.Y()-ak.X()*ai.Y())  );

		k = vec3<T>(	  ai.Y()*aj.Z()-aj.Y()*ai.Z(),
						-(ai.X()*aj.Z()-aj.X()*ai.Z()),
						  ai.X()*aj.Y()-aj.X()*ai.Y()  );

		scale(1.0/det);
		transpose();
	}
	void inverse(const mat3<T> &a)
	{
		mat3<T> r(a);
		T dd = r.det();
		if (dd==T(0)) { printf("cannot inverse with nul det!"); exit(1);  }
		inverse(r, dd);
	}
	void inverse(T dd)
	{
		vec3<T> ai(i), aj(j), ak(k);
		if (dd==T(0)) { printf("cannot inverse with nul det!"); exit(1);  }
		i = vec3<T>(	  aj.Y()*ak.Z()-ak.Y()*aj.Z(),
						-(aj.X()*ak.Z()-ak.X()*aj.Z()),
						  aj.X()*ak.Y()-ak.X()*aj.Y()  );

		j = vec3<T>(	-(ai.Y()*ak.Z()-ak.Y()*ai.Z()),
						  ai.X()*ak.Z()-ak.X()*ai.Z(),
						-(ai.X()*ak.Y()-ak.X()*ai.Y())  );

		k = vec3<T>(	  ai.Y()*aj.Z()-aj.Y()*ai.Z(),
						-(ai.X()*aj.Z()-aj.X()*ai.Z()),
						  ai.X()*aj.Y()-aj.X()*ai.Y()  );

		scale(1.0/dd);
		transpose();
	}
	void inverse()
	{
	T dd = det();
	if (dd==T(0)) { printf("cannot inverse with nul det!"); exit(1); }
	inverse(dd);
	}
	// transpose matrix (symetry along the diagonal)
	void transpose()
	{
		mat3<T> r;

		r.i=vec3<T>(i.X(), j.X(), k.X());
		r.j=vec3<T>(i.Y(), j.Y(), k.Y());
		r.k=vec3<T>(i.Z(), j.Z(), k.Z());
		*this =  r;
	}

};

#endif // !defined(AFX_MAT3_H__1544EA11_662A_41FD_8BD3_D7311F73A131__INCLUDED_)

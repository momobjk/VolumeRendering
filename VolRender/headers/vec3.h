// vec3.h: interface for the generic vec3 class.
//
// defines a vector in 3D
// By JMD 8/8/04
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VEC3_H__A6721A04_8090_4DA7_A7BB_8AAB841F9C11__INCLUDED_)
#define AFX_VEC3_H__A6721A04_8090_4DA7_A7BB_8AAB841F9C11__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <cmath>
#include <stdlib.h>
#include <cstdio>

template <class T> class vec3  
{
protected:
	T x, y, z;

public:
	//////// constructors
	vec3<T>():x(T(0)),y(T(0)),z(T(0)) {  }
	vec3<T>(T a):x(a),y(a),z(a) {  }
	vec3<T>(T a, T b, T c):x(a),y(b),z(c) { }

	// define vector from two points (a,b): v=ab
	void PVec(const vec3<T> &a, const vec3<T> &b) { x=b.x-a.x; y=b.y-a.y; z=b.z-a.z; }

	//////// selectors of coordinates
	T X() const { return x; }
	T Y() const { return y; }
	T Z() const { return z; }
	T &operator[](int i)
	{
		if (i<0 || i>2) { printf("out of vec3 range!"); }
		if (i==0) return x;
		if (i==1) return y;
		return z;
	}


	///////  operations
	void add(const vec3<T> &a, const vec3<T> &b) { x=a.x+b.x; y=a.y+b.y; z=a.z+b.z; }
	vec3<T> operator+(const vec3<T> &b) const { vec3<T> vv; vv.x=x+b.x; vv.y=y+b.y; vv.z=z+b.z; return vv;}
	void operator+=(const vec3<T> &b) { x+=b.x; y+=b.y; z+=b.z; }

	void sub(const vec3<T> &a, const vec3<T> &b) { x=a.x-b.x; y=a.y-b.y; z=a.z-b.z; }
	vec3<T> operator-(const vec3<T> &b) const { vec3<T> vv; vv.x=x-b.x; vv.y=y-b.y; vv.z=z-b.z; return vv;}
	void operator-=(const vec3<T> &b) { x-=b.x; y-=b.y; z-=b.z; }

	void mult(const vec3<T> &a, const vec3<T> &b) { x=b.x*a.x; y=b.y*a.y; z=b.z*a.z; }
	vec3<T> operator*(const vec3<T> &b) const { vec3<T> vv; vv.x=x*b.x; vv.y=y*b.y; vv.z=z*b.z; return vv;}
	void operator*=(const vec3<T> &b) { x*=b.x; y*=b.y; z*=b.z; }
	void operator*=(T v) { x*=v; y*=v; z*=v; }

	void div(const vec3<T> &a, const vec3<T> &b) { x=a.x/b.x; y=a.y/b.y; z=a.z/b.z; }
	vec3<T> operator/(const vec3<T> &b) const { vec3<T> vv; vv.x=x/b.x; vv.y=y/b.y; vv.z=z/b.z; return vv;}
	void operator/=(const vec3<T> &b) { x/=b.x; y/=b.y; z/=b.z; }
	void operator/=(T v) { x/=v; y/=v; z/=v;}

	// scale by scalar value
	void scale(const vec3<T> &a, T k) { x=a.x*k; y=a.y*k; z=a.z*k; }
	void scale(T k) { x*=k; y*=k; z*=k; }
	// divide by scalar, this scalar cannot be null
	void normalize(T norme) { x /=  norme; y /=  norme; z /=  norme; }
	void divScale(T k) { x/=k; y/=k; z/=k; }

	// absolute value component by component
	void abs()
		{
		if (x<T(0)) x = -x;
		if (y<T(0)) y = -y;
		if (z<T(0)) z = -z;
		}

	// Trunc to max value
	void trunc(T max)
		{
		if (x>max) x = max;
		if (y>max) y = max;
		if (z>max) z = max;
		}

	// random vector
	void random()
	{
		x = ((T)rand()/(T)RAND_MAX);
		y = ((T)rand()/(T)RAND_MAX);
		z = ((T)rand()/(T)RAND_MAX);
	}

	// points
	vec3<T> shift(const vec3<T> &p, const vec3<T> &dir, T length)
	{
		x = p.x+length*dir.x;
		y = p.y+length*dir.y;
		z = p.z+length*dir.z;
	}
	void boxPoints(const vec3<T> &dir, T length, vec3<T> points[]) const
	{
		vec3<T> min = *this;
		vec3<T> max = dir;
		max.scale(length);
		max.add(max, min);
		points[0]=min;
		points[1]=vec3<T>(min.X(),min.Y(),max.Z());
		points[2]=vec3<T>(min.X(),max.Y(),min.Z());
		points[3]=vec3<T>(min.X(),max.Y(),max.Z());
		points[4]=vec3<T>(max.X(),min.Y(),min.Z());
		points[5]=vec3<T>(max.X(),min.Y(),max.Z());
		points[6]=vec3<T>(max.X(),max.Y(),min.Z());
		points[7]=max;
	}


	// dot product (produit scalaire)
	T dot(const vec3<T> &v) const { return x*v.x+y*v.y+z*v.z; }

	// opposite vector
	void reverse() { x=-x; y=-y; z=-z; }
	vec3<T> operator-() { vec3<T> v(-x,-y,-z); return v; }

	// cross product (produit vectoriel)
	void cross( const vec3<T> &v1, const vec3<T> &v2)
		{
		  vec3<T> v;

		  v.x =   (v1.y)*(v2.z)  - ((v1.z)*(v2.y));
		  v.y = -((v1.x)*(v2.z)) + ((v1.z)*(v2.x));
		  v.z =   (v1.x)*(v2.y)  - ((v1.y)*(v2.x));
		  *this = v;
		 } 

	// Makes a rotation of vector v according to x, y and z axis sequentially
	//vector a indicates the respective angles in radian (0-2PI)
	void rotateCoord( const vec3<T> &v, const vec3<T> &a)
		{
		vec3<T> p,p2;

		/* Rotation along X axis */
		p.z = v.z*cos(a.x)-v.y*sin(a.x); 
		p.y = v.z*sin(a.x)+v.y*sin(a.x);
		p.x = v.x;

		/* Rotation along Y axis */
		p2.x = p.x*cos(a.y)-p.z*sin(a.y);                                
		p2.z = p.x*sin(a.y)+p.z*cos(a.y);
		p2.y = p.y;

		/* Rotation along Z axis */
		x = p2.x*cos(a.z)-p2.y*sin(a.z);                                 
		y = p2.x*sin(a.z)+p2.y*cos(a.z);
		z = p2.z;
		}

	// Makes a rotation of vector v according to x, y and z axis sequentially
	// operates like rotateCoord, but this time cosine and sine are directly given by vectors
	void rotateCosSinCoord( const vec3<T> &v, const vec3<T> &co, const vec3<T> &si)
		{
		vec3<T> p,p2;

		/* Rotation along X axis */
		p.z = v.z*co.x-v.y*si.x;
		p.y = v.z*si.x+v.y*co.x;
		p.x = v.x;

		/* Rotation along Y axis */
		p2.x = p.x*co.y-p.z*si.y;
		p2.z = p.x*si.y+p.z*co.y;
		p2.y = p.y;

		/* Rotation along Z axis */
		x = p2.x*co.z-p2.y*si.z;
		y = p2.x*si.z+p2.y*co.z;
		z = p2.z;
		}

	// chooses an orthogonal vector to r, parallel to the z=0 plane 
	//the resulting vector is normalized 
	void orthogonal(const vec3<T> &r)
		{ 
		T no;

		if (r.x==0 && r.y==0) { x=0.0; y=1.0; }
		else { x = -r.y; y = r.x; }
		z = 0.0;
		if (y<0.0) reverse();
		else if (x<0.0) reverse();
		no=norm(); if (no==T(0)) { x=T(0); y=T(0); z=T(0); }
		normalize(no);
		}
	// chooses an orthogonal vector to r, parallel to plane defined by the couple (r,v) 
	//the resulting vector is normalized 
	void orthogonalInPlan(const vec3<T> &r, const vec3<T> &v)
		{ 
		*this = r; scale(-r.dot(v)/r.dot(r));
		this->operator+=(v);
		T no=norm(); if (no==T(0)) { x=T(0); y=T(0); z=T(0); }
		normalize(no);
		}

	// compute median vector between a and b, result is not normalized
	void bissec( const vec3<T> &a, const vec3<T> &b)
		{
		add(a,b); scale(0.5);
		}

	// compute reflection vector from normal n and viewing direction v, 
	//dot product between n and v must be positive
	//Note: result is not normalized 
	void reflection( const vec3<T> &v, const vec3<T> &n)
		{
		T scal;
		vec3<T> r;

		scal= 2.0*n.dot(v);
		r.scale(n,scal); 
		r.sub(r,v);
		*this=r;
		}

	// compute refraction vector from normal n and viewing direction v
	//according to the refraction indice ir (ratio between the two media) 
	//and Snell's law,
	//dot product between n and v must be positive
	//result is normalized,
	//returns false if computation was not possible (result is complex number), in this
	//case resulting vector is set null
	bool refraction( const vec3<T> &v, const vec3<T>  &n, double ir)
		{
		vec3<T>  vp, xx, r;
		double	kf;

		if (ir==1.0) { *this=v; reverse(); kf=(double)norm(); normalise((T)kf); return true; }
		vp = v;
		kf = -1.0 / v.dot(n);  
		vp.scale(kf);
		xx.add(n, vp);
		kf = ir*ir*vp.dot(vp) - xx.dot(xx);
		if (kf <= T(0)) {  x=T(0); y=T(0); z=T(0); return false; }
		kf = 1.0 / sqrt(kf);
		r = xx; 
		r.scale(T(kf));
		r.sub(r,n); 
		kf=r.norm();
		r.normalize(kf);
		*this = r;
		return true;
		}



	// linear interpolation between v1 and v2 
	//according to t (between 0 and 1),
	//if t=0 result is v1, if t=1 result is v2
	//result is not normalized 
	void interpolate(const vec3<T> &v1, const vec3<T> &v2, T t)
		{
		x = (1.0-t)*v1.x + t*v2.x; 
		y = (1.0-t)*v1.y + t*v2.y; 
		z = (1.0-t)*v1.z + t*v2.z; 
		}

	// linear interpolation between v1 and v2 
	//according to vector t; each component of t is between 0 and 1 and gives coefficient,
	//result is not normalized 
	void interpolate(const vec3<T> &v1, const vec3<T> &v2, const vec3<T> &t)
		{
		x = (1.0-t.x)*v1.x + t.x*v2.x; 
		y = (1.0-t.y)*v1.y + t.y*v2.y; 
		z = (1.0-t.z)*v1.z + t.z*v2.z; 
		}


	// compares component by component and keeps the min one 
	//result is not normalized
	void keepMin(const vec3<T> &v1, const vec3<T> &v2)
		{
		if (v1.x<v2.x) x = v1.x; else x = v2.x;
		if (v1.y<v2.y) y = v1.y; else y = v2.y;
		if (v1.z<v2.z) z = v1.z; else z = v2.z;
		}

	// compares component by component and keeps the max one 
	//result is not normalized 
	void keepMax(const vec3<T> &v1, const vec3<T> &v2)
		{
		if (v1.x>v2.x) x = v1.x; else x = v2.x;
		if (v1.y>v2.y) y = v1.y; else y = v2.y;
		if (v1.z>v2.z) z = v1.z; else z = v2.z;
		}

	
	// returns the largest / smallest component
	T maxCoord() const
		{
		if (x>y && x>z) return x;
		if (y>z) return y;
		return z;
		}
	T minCoord() const
		{
		if (x<y && x<z) return x;
		if (y<z) return y;
		return z;
		}


	// returns the euclidian distance between two points
	T distance( const vec3<T> &a) const
		{
		vec3<T> r;

		r.PVec(a, *this);
		return r.norm(); 
		}

	// compares two vectors
	bool equals(const vec3<T> &p1) const
		{
		return (p1.x==x) && (p1.y==y) && (p1.z==z);
		}
	bool operator==(const vec3<T> &p1) const
		{
		return (p1.x==x) && (p1.y==y) && (p1.z==z);
		}
	// true if the points are closer than dist
	bool areClose(const vec3<T> &p1, T dist) const
		{
		vec3<T> v;
		v.PVec(p1, *this);
		v.abs();
		if (v.x<dist && v.y<dist && v.z<dist) return true;
		return false;
		}

	// computes the norm (vector length) 
	T normSquared() const { return x*x+y*y+z*z; } 
	T norm() const { return sqrt(normSquared()); }

	// compares to null vector
	bool isNull() const
		{
		return x==T(0) && y==T(0) && z==T(0);
		}

	bool isQNull(T dist) const
		{
		vec3<T> v=*this;
		v.abs();
		return v.x<=dist && v.y<=dist && v.z<=dist;
		}

};


#endif // !defined(AFX_VEC3_H__A6721A04_8090_4DA7_A7BB_8AAB841F9C11__INCLUDED_)

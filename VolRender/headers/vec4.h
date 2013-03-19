//vec4.h: interface for the vec4 class.
//
// Defines a vector in 4D homogenous coordinates
//
// By JMD 7/8/04
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VEC4_H__483CB07F_3D52_461D_9ADD_5849767E0AB7__INCLUDED_)
#define AFX_VEC4_H__483CB07F_3D52_461D_9ADD_5849767E0AB7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vec3.h"


template <class T> class vec4 
{
protected:
	T x, y, z, w;
public:
	//////// constructors
	vec4<T>():x(T(0)), y(T(0)), z(T(0)), w(T(1)) { }
	vec4<T>(T a):x(a),y(a),z(a),w(T(1)) {  }
	vec4<T>(T a, T b, T c):x(a),y(b),z(c),w(T(1)) { }
	vec4<T>(T a, T b, T c, T d):x(a),y(b),z(c),w(d) {  }
	vec4<T>(const vec3<T> &a):x(a.X()),y(a.Y()),z(a.Z()),w(T(1)) {  }
	vec4<T>(const vec3<T> &a, T d):x(a.X()),y(a.Y()),z(a.Z()),w(d) { }

	// cast into a 3D vector, works only if w component is not 0
	operator vec3<T>()
		{
		return vec3<T>(x/w,y/w,z/w);
		}

	//////// selectors
	T X()const { return x; }
	T Y()const { return y; }
	T Z()const { return z; }
	T W()const { return w; }

	T &operator[](int i)
	{
		if (i<0 || i>3) { printf("out of vec4 range!"); }
		if (i==0) return x;
		else if (i==1) return y;
		else if (i==2) return z;
		return w;
	}

	///////  operations
	void add(const vec4<T> &a, const vec4<T> &b) { x=a.x+b.x; y=a.y+b.y; z=a.z+b.z; w=a.w+b.w; }
	vec4<T> operator+(const vec4<T> &b) const { vec4<T> vv; vv.x=x+b.x; vv.y=y+b.y; vv.z=z+b.z; vv.w=w+b.w; return vv;}
	void operator+=(const vec4<T> &b) { x+=b.x; y+=b.y; z+=b.z; w+=b.w; }

	void sub(const vec4<T> &a, const vec4<T> &b) { x=a.x-b.x; y=a.y-b.y; z=a.z-b.z; w=a.w-b.w; }
	vec4<T> operator-(const vec4<T> &b) const { vec4<T> vv; vv.x=x-b.x; vv.y=y-b.y; vv.z=z-b.z; vv.w=w-b.w; return vv;}
	void operator-=(const vec4<T> &b) { x-=b.x; y-=b.y; z-=b.z; w-=b.w; }

	void mult(const vec4<T> &a, const vec4<T> &b) { x=b.x*a.x; y=b.y*a.y; z=b.z*a.z; w=b.w*a.w; }
	vec4<T> operator*(const vec4<T> &b) const { vec4<T> vv; vv.x=x*b.x; vv.y=y*b.y; vv.z=z*b.z; vv.w=w*b.w; return vv;}
	void operator*=(const vec4<T> &b) { x*=b.x; y*=b.y; z*=b.z; w*=b.w; }
	void operator*=(T v) { x*=v; y*=v; z*=v; w*=v; }

	void div(const vec4<T> &a, const vec4<T> &b) { x=a.x/b.x; y=a.y/b.y; z=a.z/b.z; w=a.w/b.w; }
	vec4<T> operator/(const vec4<T> &b) const { vec4<T> vv; vv.x=x/b.x; vv.y=y/b.y; vv.z=z/b.z; vv.w=w/b.w; return vv;}
	void operator/=(const vec4<T> &b) { x/=b.x; y/=b.y; z/=b.z; w/=b.w; }
	void operator/=(T v) { x/=v; y/=v; z/=v; w/=v; }

	// scale by scalar value
	void scale(const vec4<T> &a, T k) { x=a.x*k; y=a.y*k; z=a.z*k; w=a.w*k; }
	void scale(T k) { x*=k; y*=k; z*=k; w*=k; }

	// absolute value component by component
	void abs()
		{
		if (x<0) x = -x;
		if (y<0) y = -y;
		if (z<0) z = -z;
		if (w<0) w = -w;
		}

	// dot product (produit scalaire)
	T dot(const vec4<T> &v) const { return x*v.x+y*v.y+z*v.z+w*v.w; }

	// opposite vector
	void reverse() { x=-x; y=-y; z=-z; w=-w; }
	vec4<T> operator-() { vec4<T> v(-x,-y,-z,-w); return v; }


	// linear interpolation between v1 and v2 
	//according to t (between 0 and 1),
	//if t=0 result is v1, if t=1 result is v2
	//result is not normalized 
	void interpolate(const vec4<T> &v1, const vec4<T> &v2, T t)
		{
		x = (1.0-t)*v1.x + t*v2.x; 
		y = (1.0-t)*v1.y + t*v2.y; 
		z = (1.0-t)*v1.z + t*v2.z; 
		w = (1.0-t)*v1.w + t*v2.w; 
		}

	// linear interpolation between v1 and v2 
	//according to vector t; each component of t is between 0 and 1 and gives coefficient,
	//result is not normalized 
	void interpolate(const vec4<T> &v1, const vec4<T> &v2, const vec4<T> &t)
		{
		x = (1.0-t.x)*v1.x + t.x*v2.x; 
		y = (1.0-t.y)*v1.y + t.y*v2.y; 
		z = (1.0-t.z)*v1.z + t.z*v2.z; 
		w = (1.0-t.w)*v1.w + t.w*v2.w; 
		}


	/* compares two vectors */
	bool equals(const vec4<T> &p1) const
		{
		return (p1.x==x) && (p1.y==y) && (p1.z==z) && (p1.w==w);
		}
	bool operator==(const vec4<T> &p1) const
		{
		return (p1.x==x) && (p1.y==y) && (p1.z==z) && (p1.w==w);
		}

};

#endif // !defined(AFX_VEC4_H__483CB07F_3D52_461D_9ADD_5849767E0AB7__INCLUDED_)

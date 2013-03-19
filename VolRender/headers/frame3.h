// frame3.h: interface for the frame3 class.
//
// frame3: defines a frame for an euclidian 3D space
// a frame is composed of an origin point and 3 orthonormal vectors
// frame3 is composed of a vec3<T> and a mat3<T>
// main operations allow to change points and vectors from one frame to another
// by extracting transfer matrices
//
// By JMD 10/8/06
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FRAME3_H__824D3162_521A_4FD7_9D86_19B3BA026926__INCLUDED_)
#define AFX_FRAME3_H__824D3162_521A_4FD7_9D86_19B3BA026926__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "mat4.h"


template <class T> class frame3 
{
protected:

	vec3<T> or;
	mat3<T> mat;

public:
	// constructors
	frame3<T>():mat(), or(vec3<T>(T(0),T(0),T(0))) { }
	frame3(const vec3<T> &o, const vec3<T> &x, const vec3<T> &y, const vec3<T> &z) : mat(x,y,z),or(o) { }
	frame3(const vec3<T> &o):mat(), or(o) { }
	frame3(const vec3<T> &o, const mat3<T> &m):mat(m), or(o) { }

	// selectors
	vec3<T> origin() const { return or; }
	vec3<T> I() const { return mat.I(); }
	vec3<T> J() const { return mat.J(); }
	vec3<T> K() const { return mat.K(); }

	// computes coordinates X,Y,Z of point pt in the frame
	vec3<T> coordinates(const vec3<T> &pt) const
	{
		vec3<T> p; p.PVec(or,pt);
		return vec3<T>(p.dot(I())/I().normSquared(),p.dot(J())/J().normSquared(),p.dot(K())/K().normSquared());
	}

	T det() const { return mat.det(); }
	bool operator==(const frame3<T> &f) const
	{
		return or==f.or && mat==f.mat;
	}
	// Compute change of frame matrix (dimension 3) for vectors
	// Resulting matrix allows to express a vector given in base coordinates :
	// (1,0,0);(0,1,0);(0,0,1) into the frame coordinates  
	mat3<T> changeToFrame3(T det) const
	{
		hvMat3<T> r;
		r.inverse(mat, det);
		return r;
	}

	// Compute change of frame matrix (dimension 3) for vectors
	// Resulting matrix allows to express a vector given in frame coordinates
	// into basic coordinates, e.g. (1,0,0);(0,1,0);(0,0,1)   
	mat3<T> changeFromFrame3() const 
	{
		return mat;
	}

	// Compute change of frame matrix (dimension 4) for points
	// Resulting matrix allows to express a point given in base coordinates :
	// O(0,0,0);i(1,0,0);j(0,1,0);k(0,0,1) into the frame coordinates  
	mat4<T> changeToFrame4(T det) const 
	{
		vec3<T> io=or; io.reverse();
		mat4<T> m1(mat3<T>(), io);
		mat3<T> r=changeToFrame3(det);
		mat4<T> m2(r);
		m1.mult(m1, m2);
		return m1;
	}

	// Compute change of frame matrix (dimension 4) for points
	// Resulting matrix allows to express a point given in frame coordinates
	// into basic coordinates, e.g. O(0,0,0);i(1,0,0);j(0,1,0);k(0,0,1)   
	mat4<T> changeFromFrame4() const
	{
		mat3<T> r=changeFromFrame3();
		mat4<T> m1(r);
		vec3<T> io=or;
		mat4<T> m2(mat3<T>(), io);
		m1.mult(m1, m2);
		return m1;
	}

};


#endif // !defined(AFX_FRAME3_H__824D3162_521A_4FD7_9D86_19B3BA026926__INCLUDED_)

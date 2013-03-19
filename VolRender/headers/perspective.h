// perspect.h: interface for the perspective class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PERSPECT_H__BB034B5A_4CC1_4FBC_B131_9AECA8B630DB__INCLUDED_)
#define AFX_PERSPECT_H__BB034B5A_4CC1_4FBC_B131_9AECA8B630DB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "frame3.h"
#include "transform.h"
#include "line3.h"

// Defines a perspective projection in 3D space: a linear transform that changes a 3D point into another 3D point belonging to a projection
// plan. The result of projection is a point P(x,y,z) where x,y,z are normalized between -1 and 1
// (x,y) defines a 2D vector according to a 2D frame located on the projection plane
// z is a distance measurement of the projected point towards the projection plane (it is normalized according
// to a maximal distance called farplane)

template <class T> class perspective : public transform<T> 
{
protected:
	frame3<T>			eyeframe;
	vec3<T>				lookat;
	T					width, height, focus, farplane;

public:
	perspective<T>():eyeframe(),lookat() { width=T(0.0); height=T(0.0); focus=T(0.0); farplane=T(1.0); }

	perspective<T>(vec3<T> &eye, const vec3<T> &at, const vec3<T> &up, T w, T h, T n, T f)
	{
		vec3<T> zaxis; zaxis.PVec(eye, at);
		zaxis.normalize(zaxis.norm());
		vec3<T> xaxis;
		xaxis.cross(up, zaxis);
		xaxis.normalize(xaxis.norm());
		zaxis.reverse();
		eyeframe = frame3<T>(eye, xaxis, up, zaxis);
		lookat = at;
		width = w; height = h; focus = n; farplane = f;
		frustum(-w/T(2.0), w/T(2.0), -h/T(2.0), h/T(2.0), n, f);
		lookAt(eye, at, up);
	}

	// selectors
	vec3<T> getOrigin() const { return eyeframe.origin(); }
	frame3<T> getFrame() const { return eyeframe; }
	vec3<T> getLookAt() const { return lookat; }

	bool project(const vec3<T> &p, vec3<T> &res) const
	{
		vec3<T> pv; pv.PVec(eyeframe.origin(), p);
		if (pv.dot(eyeframe.K())>=0.0) return false;
		vec4<T> u(p);
		vec4<T> v = mult(u);
		if (v.W()==T(0.0)) return false;
		res = vec3<T>(v);
		return true;
	}
	
	line3<T> reverse(const vec3<T> &p) const
	{
		vec3<T> r=eyeframe.K();
		r.scale(-focus);
		vec3<T> xx=eyeframe.I();
		xx.scale(p.X()*width/T(2.0));
		vec3<T> yy=eyeframe.J();
		yy.scale(p.Y()*height/T(2.0));
		r+=xx; r+=yy;
		r.normalize(r.norm());
		return line3<T>(eyeframe.origin(), r);
	}
	
};


#endif // !defined(AFX_MONOPERSPECT_H__BB034B5A_4CC1_4FBC_B131_9AECA8B630DB__INCLUDED_)

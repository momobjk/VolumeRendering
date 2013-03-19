// line3.h: interface for the line3 class.
//
// line3 defines a line in a 3D space: L= O+Dl, where O is origin and D direction
// 
// By JMD 10/8/06
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LINE3_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_)
#define AFX_LINE3_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vec3.h"

template <class T> class line3
{
protected:
	vec3<T>		o;	// origin
	vec3<T>		d;	// direction
public:
	// constructors
	line3<T>() { }

	// defines a line by a point pt (origin) and a vector dir (direction)
	// dir must be normalized
	line3<T>(const vec3<T> &pt, const vec3<T> &dir) { o=pt; d=dir; }

	// defines a line from origin (0,0,0) in direction dir
	// dir must be normalized
	line3<T>(const vec3<T> &dir):o(T(0)) { d=dir; }

	// selectors
	vec3<T> origin() const { return o; }
	vec3<T> direction() const { return d; }

	// Compute a point (result) on the line at distance l from the origin
	vec3<T> pointOfLine(T l) const
	{
		vec3<T> p; p.scale(d, l);
		p += o;
		return p;
	}

	void reverse() { d.reverse(); }

};

#endif // !defined(AFX_LINE3_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_)

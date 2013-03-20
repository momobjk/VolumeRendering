#ifndef RAW3DDATA_H
#define RAW3DDATA_H

#include "vec3.h"
#include <cstdio>

class Raw3DData
{
    private:
        int             sx,sy,sz;	// size of data array
        unsigned char   *data;   // raw data array
        unsigned char   *grad;   // gradient data for lighting on GPU - to be computed or loaded
        double          scale;  // for normalization
        unsigned char	*transfer;
        unsigned short	*preIntTab,*rampLeft,*rampRight;  // pre-integration tables
        int             tabSize;
        double			voxelStep;

    public :
        Raw3DData();
        unsigned char get(int ind[3]) const;
        unsigned char get(int i, int j, int k) const;
        double rawdensity(int i, int j, int k, double ri, double rj, double rk) const;
        vec3<double> rawgradient(int i, int j, int k, double ri, double rj, double rk) const;
        // Loading RAW 3D texture images on 8bits
        bool importRAW(const char *name, int xx, int yy, int zz, double ss);
        bool importGrad(const char *name);
        void computeGrad();
        void computePreIntegration(unsigned char *tf, int ts, double vs, int nstep);

        //getters
        int getSx() const{ return sx; }
        int getSy() const{ return sy; }
        int getSz() const{ return sz; }
        unsigned char* getData() { return data; }
        unsigned char* getGrad() { return grad; }
        double getScale() { return scale; }
        unsigned char* getTransfer() { return transfer; }
        unsigned short* getPreIntTab() { return preIntTab; }
        unsigned short* getRampLeft() { return rampLeft; }
        unsigned short* getRampRight() { return rampRight; }
        int getTabSize() { return tabSize; }
        double getVoxelStep() { return voxelStep; }
};


#endif // RAW3DDATA_H

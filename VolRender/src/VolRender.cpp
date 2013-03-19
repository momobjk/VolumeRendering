// VolRender.cpp : Defines the entry point for the console application.
//

#include <cstdio>
#include <cmath>
#include <cstring>
#include <GL/glew.h>
#include <GL/glut.h>

#include "transform.h"

#define PI 3.1415926535897932384626433832795

//#define PREINT

// GL ERROR CHECK
int CheckGLError(char *file, int line)
{
    GLenum glErr;
    int    retCode = 0;

    glErr = glGetError();
    while (glErr != GL_NO_ERROR) {
        printf("GL Error #%d ( %s ) in File %s at line: %d\n", glErr, gluErrorString(glErr),file, line );
        retCode = 1;
        glErr = glGetError();
    }
    return retCode;
}
#define CHECK_GL_ERROR() CheckGLError(__FILE__, __LINE__)

class Raw3DData 
{
public:
    int					sx,sy,sz;	// size of data array
    unsigned char		*data;   // raw data array
    unsigned char		*grad;   // gradient data for lighting on GPU - to be computed or loaded
    double				scale;  // for normalization

    unsigned char	*transfer;
    unsigned short	*preinttab,*rampleft,*rampright;  // pre-integration tables
    int					tabsize;
    double			voxelstep;

    Raw3DData() { sx=0; sy=0; sz=0; data=0; grad=0; scale = 1.0/255.0; }
    unsigned char get(int ind[3]) const { return data[ind[0]+ind[1]*sx+ind[2]*sx*sy]; }
    unsigned char get(int i, int j, int k) const { return data[i+j*sx+k*sx*sy]; }
    double rawdensity(int i, int j, int k, double ri, double rj, double rk) const
    {
        if (i>=sx-1 || j>=sy-1 || k>=sz-1 || i<0 || j<0 || k<0) { return 0.0; }
        int id1[3]; id1[0]=i; id1[1]=j; id1[2]=k;
        int id2[3]; id2[0]=i+1; id2[1]=j; id2[2]=k;
        double vj1 = (double)get(id1)*scale*(1.0-(double)ri)+(double)get(id2)*scale*(double)ri;
        id1[0]=i; id1[1]=j+1; id1[2]=k;
        id2[0]=i+1; id2[1]=j+1; id2[2]=k;
        double vj2 = (double)get(id1)*scale*(1.0-(double)ri)+(double)get(id2)*scale*(double)ri;
        id1[0]=i; id1[1]=j; id1[2]=k+1;
        id2[0]=i+1; id2[1]=j; id2[2]=k+1;
        double vk1 = (double)get(id1)*scale*(1.0-(double)ri)+(double)get(id2)*scale*(double)ri;
        id1[0]=i; id1[1]=j+1; id1[2]=k+1;
        id2[0]=i+1; id2[1]=j+1; id2[2]=k+1;
        double vk2 = (double)get(id1)*scale*(1.0-(double)ri)+(double)get(id2)*scale*(double)ri;
        vj1 = vj1*(1.0-(double)rj)+vj2*(double)rj;
        vk1 = vk1*(1.0-(double)rj)+vk2*(double)rj;
        return double(vj1*(1.0-(double)rk)+vk1*(double)rk);
    }
    vec3<double> rawgradient(int i, int j, int k, double ri, double rj, double rk) const
    {
        double dx = rawdensity(i-1, j, k, ri, rj, rk);
        dx -= rawdensity(i+1, j, k, ri, rj, rk);
        double dy = rawdensity(i, j-1, k, ri, rj, rk);
        dy -= rawdensity(i, j+1, k, ri, rj, rk);
        double dz = rawdensity(i, j, k-1, ri, rj, rk);
        dz -= rawdensity(i, j, k+1, ri, rj, rk);
        vec3<double> vgrad = vec3<double>(dx,dy,dz);
        if (vgrad.norm()==double(0)) return vec3<double>();
        vgrad.normalize(vgrad.norm());
        return vgrad;
    }
    // Loading RAW 3D texture images on 8bits
    bool importRAW(char *name, int xx, int yy, int zz, double ss)
    {

        FILE *fd = fopen(name, "rb");
        if (fd==0) { printf("cannot open file!\n"); exit(1); }
        sx=xx; sy=yy; sz=zz;
        data = new unsigned char [sx*sy*sz];
        int i,j,k;
        unsigned char *dd = new unsigned char [sx];
        for (k=0; k<sz; k++)
            for (j=0; j<sy; j++)
            {
                fread(dd, sizeof(unsigned char), sx, fd);
                for (i=0; i<sx; i++)
                {
                    data[i+j*sx+k*sx*sy] = dd[i];
                }
            }
        scale = ss;
        delete [] dd;
        return true;
    }
    bool importGrad(char *name)
    {

        FILE *fd = fopen(name, "rb");
        if (fd==0) { printf("cannot open file!\n"); exit(1); }
        grad = new unsigned char [4*sx*sy*sz];
        fread(grad,sx*sy*sz*4, sizeof(unsigned char), fd);
        fclose(fd);
        return true;
    }
    void computeGrad()
    {
        int i,j,k;
        if (data==0) return;
        grad = new unsigned char [4*sx*sy*sz];
        printf("computing gradient.\n");
        for (i=0; i<sx; i++) for (j=0; j<sy; j++) for (k=0; k<sz; k++)
        {
            //printf("(%d, %d, %d)",i,j,k);
            vec3<double> vgrad = rawgradient(i,j,k,0.0,0.0,0.0);
            double dens = rawdensity(i,j,k,0.0,0.0,0.0);
            int ind = i+j*sx+k*sx*sy;
            grad[4*ind]=(unsigned char)(127.0*vgrad.X()+127.0);
            grad[4*ind+1]=(unsigned char)(127.0*vgrad.Y()+127.0);
            grad[4*ind+2]=(unsigned char)(127.0*vgrad.Z()+127.0);
            grad[4*ind+3]=(unsigned char)(255.0*dens);
        }
        printf("done!\n");
    }
    void computePreIntegration(unsigned char *tf, int ts, double vs, int nstep)
    {
        int i,j,k;
        double scal=65535.0;

        tabsize=ts;
        voxelstep=vs;
        transfer = tf;

        preinttab=(unsigned short *)malloc(ts*ts*4*sizeof (unsigned short));
        rampleft=(unsigned short *)malloc(ts*ts*4*sizeof (unsigned short));
        rampright=(unsigned short *)malloc(ts*ts*4*sizeof (unsigned short));
        for (i=0; i<tabsize*tabsize; i++) { preinttab[i]=0; rampleft[i]=0; rampright[0]=0; }

        for (i=0; i<tabsize; i++)
        {
            double alpha=0.0, colt=1.0, alphamid=0.0;
            double vr = 0.0, vg = 0.0, vb = 0.0, alpha0r=0.0;
            double vrl = 0.0, vgl = 0.0, vbl = 0.0, alphal=0.0;
            double vrr = 0.0, vgr = 0.0, vbr = 0.0, alphar=0.0;
            double vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
            double vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;

            double vrg = 0.0, vgg = 0.0, vbg = 0.0;
            double vrlg = 0.0, vglg = 0.0, vblg = 0.0;
            double vr0rg = 0.0, vg0rg = 0.0, vb0rg = 0.0;
            double vrrg = 0.0, vgrg = 0.0, vbrg = 0.0;
            double vr2rg = 0.0, vg2rg = 0.0, vb2rg = 0.0;
            double vr3rg = 0.0, vg3rg = 0.0, vb3rg = 0.0;

            double colr, colg,colb, de;
            //printf("step %d/%d\n", i,pi.sizeX());
            alpha=0.0; colt=1.0; alphamid=0.0; vr = 0.0; vg = 0.0; vb = 0.0; alpha0r=0.0;
            vrl = 0.0; vgl = 0.0; vbl = 0.0; alphal=0.0;
            vrr = 0.0; vgr = 0.0; vbr = 0.0; alphar=0.0;
            vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
            vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;
            vrg = 0.0; vgg = 0.0; vbg = 0.0;
            vrlg = 0.0; vglg = 0.0; vblg = 0.0;
            vrrg = 0.0; vgrg = 0.0; vbrg = 0.0;
            vr0rg = 0.0, vg0rg = 0.0, vb0rg = 0.0;
            vr2rg = 0.0, vg2rg = 0.0, vb2rg = 0.0;
            vr3rg = 0.0, vg3rg = 0.0, vb3rg = 0.0;
            colr = (double)tf[4*i+0]/255.0;
            colg = (double)tf[4*i+1]/255.0;
            colb = (double)tf[4*i+2]/255.0;
            de = (double)tf[4*i+3]/255.0;
            for (k=0; k<nstep; k++)
            {
                double tt = (double)k/(double)(nstep-1);
                alpha=colt*de/(double)nstep*voxelstep;
                vr += alpha*colr;
                vg += alpha*colg;
                vb += alpha*colb;
                alpha0r += alpha*1.0;
                vrl += alpha*colr*(1.0-tt);
                vgl += alpha*colg*(1.0-tt);
                vbl += alpha*colb*(1.0-tt);
                alphal += alpha*(1.0-tt);
                vrr += alpha*colr*tt;
                vgr += alpha*colg*tt;
                vbr += alpha*colb*tt;
                alphar+= alpha*tt;
                vr2r += alpha*colr*tt*tt;
                vg2r += alpha*colg*tt*tt;
                vb2r += alpha*colb*tt*tt;
                alpha2r += alpha*tt*tt;
                vr3r += alpha*colr*tt*tt*tt;
                vg3r += alpha*colg*tt*tt*tt;
                vb3r += alpha*colb*tt*tt*tt;
                alpha3r += alpha*tt*tt*tt;
                alphamid += alpha*(tt<0.5?tt*2.0:2.0-2.0*tt);
                colt *= (1.0-de/(double)nstep*voxelstep);
            }
            preinttab[4*(i*ts+i)]=(unsigned short)(scal*vr); preinttab[4*(i*ts+i)+1]=(unsigned short)(scal*vg); preinttab[4*(i*ts+i)+2]=(unsigned short)(scal*vb); preinttab[4*(i*ts+i)+3]=(unsigned short)(scal*(1.0-colt));
            rampleft[4*(i*ts+i)]=(unsigned short)(scal*vrl); rampleft[4*(i*ts+i)+1]=(unsigned short)(scal*vgl); rampleft[4*(i*ts+i)+2]=(unsigned short)(scal*vbl); rampleft[4*(i*ts+i)+3]=(unsigned short)(scal*alphal);
            rampright[4*(i*ts+i)]=(unsigned short)(scal*vrr); rampright[4*(i*ts+i)+1]=(unsigned short)(scal*vgr); rampright[4*(i*ts+i)+2]=(unsigned short)(scal*vbr); rampright[4*(i*ts+i)+3]=(unsigned short)(scal*alphar);
            printf("i=%d  Ka=%g,%g,%g   Kd=%g,%g,%g\n",i,vr,vg,vb,vrr,vgr,vbr);
            for (j=i+1; j<ts; j++)
            {
                double d = 1.0/(double)(j-i);
                double kk=0.0;
                alpha=0.0; colt=1.0; alphamid=0.0;  vr = 0.0; vg = 0.0; vb = 0.0; alpha0r=0.0;
                vrl = 0.0; vgl = 0.0; vbl = 0.0; alphal=0.0;
                vrr = 0.0; vgr = 0.0; vbr = 0.0; alphar=0.0;
                vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
                vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;
                vrg = 0.0, vgg = 0.0, vbg = 0.0;
                vr0rg = 0.0, vg0rg = 0.0, vb0rg = 0.0;
                vrlg = 0.0; vglg = 0.0; vblg = 0.0;
                vrrg = 0.0; vgrg = 0.0; vbrg = 0.0;
                vr2rg = 0.0, vg2rg = 0.0, vb2rg = 0.0;
                vr3rg = 0.0, vg3rg = 0.0, vb3rg = 0.0;
                for (k=0; k<nstep*(j-i); k++)
                {
                    double ramp = (double)k/(double)(nstep*(j-i)-1);
                    int ind=k/nstep;
                    double t = (double)(k-ind*nstep)/(double)nstep;
                    ind+=i;
                    colr = (double)tf[4*ind]/255.0*(1.0-t)+(double)tf[4*(ind+1)]/255.0*t;
                    colg = (double)tf[4*ind+1]/255.0*(1.0-t)+(double)tf[4*(ind+1)+1]/255.0*t;
                    colb = (double)tf[4*ind+2]/255.0*(1.0-t)+(double)tf[4*(ind+1)+2]/255.0*t;
                    de = (double)tf[4*ind+3]/255.0*(1.0-t)+(double)tf[4*(ind+1)+3]/255.0*t;
                    ind++;

                    alpha = colt*de/(double)nstep*d*voxelstep;
                    vr += alpha*colr;
                    vg += alpha*colg;
                    vb += alpha*colb;
                    alpha0r += alpha*1.0;
                    vrl += alpha*colr*(1.0-ramp);
                    vgl += alpha*colg*(1.0-ramp);
                    vbl += alpha*colb*(1.0-ramp);
                    alphal += alpha*(1.0-ramp);
                    vrr += alpha*colr*ramp;
                    vgr += alpha*colg*ramp;
                    vbr += alpha*colb*ramp;
                    alphar += alpha*ramp;
                    vr2r += alpha*colr*ramp*ramp;
                    vg2r += alpha*colg*ramp*ramp;
                    vb2r += alpha*colb*ramp*ramp;
                    alpha2r += alpha*ramp*ramp;
                    vr3r += alpha*colr*ramp*ramp*ramp;
                    vg3r += alpha*colg*ramp*ramp*ramp;
                    vb3r += alpha*colb*ramp*ramp*ramp;
                    alpha3r += alpha*ramp*ramp*ramp;
                    alphamid += alpha*(ramp<0.5?ramp*2.0:2.0-2.0*ramp);
                    colt *= (1.0-de/(double)nstep*d*voxelstep);

                }
                if (vr>1.0) vr=1.0; else if (vr<0.0) vr=0.0;
                if (vg>1.0) vg=1.0; else if (vg<0.0) vg=0.0;
                if (vb>1.0) vb=1.0; else if (vb<0.0) vb=0.0;
                if (alpha0r>1.0) alpha0r=1.0; else if (alpha0r<0.0) alpha0r=0.0;
                if (alpha>1.0) alpha=1.0; else if (alpha<0.0) alpha=0.0;
                if (vrl>1.0) vrl=1.0; else if (vrl<0.0) vrl=0.0;
                if (vgl>1.0) vgl=1.0; else if (vgl<0.0) vgl=0.0;
                if (vbl>1.0) vbl=1.0; else if (vbl<0.0) vbl=0.0;
                if (alphal>1.0) alphal=1.0; else if (alphal<0.0) alphal=0.0;
                if (vrr>1.0) vrr=1.0; else if (vrr<0.0) vrr=0.0;
                if (vgr>1.0) vgr=1.0; else if (vgr<0.0) vgr=0.0;
                if (vbr>1.0) vbr=1.0; else if (vbr<0.0) vbr=0.0;
                if (alphar>1.0) alphar=1.0; else if (alphar<0.0) alphar=0.0;
                if (vr2r>1.0) vr2r=1.0; else if (vr2r<0.0) vr2r=0.0;
                if (vg2r>1.0) vg2r=1.0; else if (vg2r<0.0) vg2r=0.0;
                if (vb2r>1.0) vb2r=1.0; else if (vb2r<0.0) vb2r=0.0;
                if (alpha2r>1.0) alpha2r=1.0; else if (alpha2r<0.0) alpha2r=0.0;
                if (vr3r>1.0) vr3r=1.0; else if (vr3r<0.0) vr3r=0.0;
                if (vg3r>1.0) vg3r=1.0; else if (vg3r<0.0) vg3r=0.0;
                if (vb3r>1.0) vb3r=1.0; else if (vb3r<0.0) vb3r=0.0;
                if (alpha3r>1.0) alpha3r=1.0; else if (alpha3r<0.0) alpha3r=0.0;
                preinttab[4*(j*ts+i)]=(unsigned short)(scal*vr); preinttab[4*(j*ts+i)+1]=(unsigned short)(scal*vg); preinttab[4*(j*ts+i)+2]=(unsigned short)(scal*vb); preinttab[4*(j*ts+i)+3]=(unsigned short)(scal*(1.0-colt));
                rampleft[4*(j*ts+i)]=(unsigned short)(scal*vrl); rampleft[4*(j*ts+i)+1]=(unsigned short)(scal*vgl); rampleft[4*(j*ts+i)+2]=(unsigned short)(scal*vbl); rampleft[4*(j*ts+i)+3]=(unsigned short)(scal*alphal);
                rampright[4*(j*ts+i)]=(unsigned short)(scal*vrr); rampright[4*(j*ts+i)+1]=(unsigned short)(scal*vgr); rampright[4*(j*ts+i)+2]=(unsigned short)(scal*vbr); rampright[4*(j*ts+i)+3]=(unsigned short)(scal*alphar);
            }
            for (j=0; j<i; j++)
            {
                double d = 1.0/(double)(i-j);
                double kk=0.0;
                colt=1.0; alpha=0.0; alphamid=0.0; vr = 0.0; vg = 0.0; vb = 0.0; alpha0r=0.0;
                vrl = 0.0; vgl = 0.0; vbl = 0.0; alphal=0.0;
                vrr = 0.0; vgr = 0.0; vbr = 0.0; alphar=0.0;
                vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
                vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;
                vrg = 0.0, vgg = 0.0, vbg = 0.0;
                vr0rg = 0.0, vg0rg = 0.0, vb0rg = 0.0;
                vrlg = 0.0; vglg = 0.0; vblg = 0.0;
                vrrg = 0.0; vgrg = 0.0; vbrg = 0.0;
                vr2rg = 0.0, vg2rg = 0.0, vb2rg = 0.0;
                vr3rg = 0.0, vg3rg = 0.0, vb3rg = 0.0;
                for (k=0; k<nstep*(i-j); k++)
                {
                    double ramp = (double)k/(double)(nstep*(i-j)-1);
                    int ind=k/nstep;
                    double t = (double)(k-ind*nstep)/(double)nstep;
                    ind+=j;
                    colr = (double)tf[4*ind]/255.0*(1.0-t)+(double)tf[4*(ind+1)]/255.0*t;
                    colg = (double)tf[4*ind+1]/255.0*(1.0-t)+(double)tf[4*(ind+1)+1]/255.0*t;
                    colb = (double)tf[4*ind+2]/255.0*(1.0-t)+(double)tf[4*(ind+1)+2]/255.0*t;
                    de = (double)tf[4*ind+3]/255.0*(1.0-t)+(double)tf[4*(ind+1)+3]/255.0*t;
                    ind++;

                    alpha = colt*de/(double)nstep*d*voxelstep;
                    vr += alpha*colr;
                    vg += alpha*colg;
                    vb += alpha*colb;
                    alpha0r += alpha*1.0;
                    vrl += alpha*colr*(1.0-ramp);
                    vgl += alpha*colg*(1.0-ramp);
                    vbl += alpha*colb*(1.0-ramp);
                    alphal += alpha*(1.0-ramp);
                    vrr += alpha*colr*ramp;
                    vgr += alpha*colg*ramp;
                    vbr += alpha*colb*ramp;
                    alphar += alpha*ramp;
                    vr2r += alpha*colr*ramp*ramp;
                    vg2r += alpha*colg*ramp*ramp;
                    vb2r += alpha*colb*ramp*ramp;
                    alpha2r += alpha*ramp*ramp;
                    vr3r += alpha*colr*ramp*ramp*ramp;
                    vg3r += alpha*colg*ramp*ramp*ramp;
                    vb3r += alpha*colb*ramp*ramp*ramp;
                    alpha3r += alpha*ramp*ramp*ramp;
                    alphamid += alpha*(ramp<0.5?ramp*2.0:2.0-2.0*ramp);
                    colt *= (1.0-de/(double)nstep*d*voxelstep);
                }
                if (vr>1.0) vr=1.0; else if (vr<0.0) vr=0.0;
                if (vg>1.0) vg=1.0; else if (vg<0.0) vg=0.0;
                if (vb>1.0) vb=1.0; else if (vb<0.0) vb=0.0;
                if (alpha0r>1.0) alpha0r=1.0; else if (alpha0r<0.0) alpha0r=0.0;
                if (alpha>1.0) alpha=1.0; else if (alpha<0.0) alpha=0.0;
                if (vrl>1.0) vrl=1.0; else if (vrl<0.0) vrl=0.0;
                if (vgl>1.0) vgl=1.0; else if (vgl<0.0) vgl=0.0;
                if (vbl>1.0) vbl=1.0; else if (vbl<0.0) vbl=0.0;
                if (alphal>1.0) alphal=1.0; else if (alphal<0.0) alphal=0.0;
                if (vrr>1.0) vrr=1.0; else if (vrr<0.0) vrr=0.0;
                if (vgr>1.0) vgr=1.0; else if (vgr<0.0) vgr=0.0;
                if (vbr>1.0) vbr=1.0; else if (vbr<0.0) vbr=0.0;
                if (alphar>1.0) alphar=1.0; else if (alphar<0.0) alphar=0.0;
                if (vr2r>1.0) vr2r=1.0; else if (vr2r<0.0) vr2r=0.0;
                if (vg2r>1.0) vg2r=1.0; else if (vg2r<0.0) vg2r=0.0;
                if (vb2r>1.0) vb2r=1.0; else if (vb2r<0.0) vb2r=0.0;
                if (alpha2r>1.0) alpha2r=1.0; else if (alpha2r<0.0) alpha2r=0.0;
                if (vr3r>1.0) vr3r=1.0; else if (vr3r<0.0) vr3r=0.0;
                if (vg3r>1.0) vg3r=1.0; else if (vg3r<0.0) vg3r=0.0;
                if (vb3r>1.0) vb3r=1.0; else if (vb3r<0.0) vb3r=0.0;
                if (alpha3r>1.0) alpha3r=1.0; else if (alpha3r<0.0) alpha3r=0.0;
                preinttab[4*(j*ts+i)]=(unsigned short)(scal*vr); preinttab[4*(j*ts+i)+1]=(unsigned short)(scal*vg); preinttab[4*(j*ts+i)+2]=(unsigned short)(scal*vb); preinttab[4*(j*ts+i)+3]=(unsigned short)(scal*(1.0-colt));
                rampleft[4*(j*ts+i)]=(unsigned short)(scal*vrl); rampleft[4*(j*ts+i)+1]=(unsigned short)(scal*vgl); rampleft[4*(j*ts+i)+2]=(unsigned short)(scal*vbl); rampleft[4*(j*ts+i)+3]=(unsigned short)(scal*alphal);
                rampright[4*(j*ts+i)]=(unsigned short)(scal*vrr); rampright[4*(j*ts+i)+1]=(unsigned short)(scal*vgr); rampright[4*(j*ts+i)+2]=(unsigned short)(scal*vbr); rampright[4*(j*ts+i)+3]=(unsigned short)(scal*alphar);
            }
        }
        for (i=0; i<ts; i++)
        {
            preinttab[4*(i*ts+0)]=0; preinttab[4*(i*ts+0)+1]=0; preinttab[4*(i*ts+0)+2]=0; preinttab[4*(i*ts+0)+3]=0;
            rampleft[4*(i*ts+0)]=0; rampleft[4*(i*ts+0)+1]=0; rampleft[4*(i*ts+0)+2]=0; rampleft[4*(i*ts+0)+3]=0;
            rampright[4*(i*ts+0)]=0; rampright[4*(i*ts+0)+1]=0; rampright[4*(i*ts+0)+2]=0; rampright[4*(i*ts+0)+3]=0;

            preinttab[4*(0*ts+i)]=0; preinttab[4*(0*ts+i)+1]=0; preinttab[4*(0*ts+i)+2]=0; preinttab[4*(0*ts+i)+3]=0;
            rampleft[4*(0*ts+i)]=0; rampleft[4*(0*ts+i)+1]=0; rampleft[4*(0*ts+i)+2]=0; rampleft[4*(0*ts+i)+3]=0;
            rampright[4*(0*ts+i)]=0; rampright[4*(0*ts+i)+1]=0; rampright[4*(0*ts+i)+2]=0; rampright[4*(0*ts+i)+3]=0;

            preinttab[4*(i*ts+ts-1)]=0; preinttab[4*(i*ts+ts-1)+1]=0; preinttab[4*(i*ts+ts-1)+2]=0; preinttab[4*(i*ts+ts-1)+3]=0;
            rampleft[4*(i*ts+ts-1)]=0; rampleft[4*(i*ts+ts-1)+1]=0; rampleft[4*(i*ts+ts-1)+2]=0; rampleft[4*(i*ts+ts-1)+3]=0;
            rampright[4*(i*ts+ts-1)]=0; rampright[4*(i*ts+ts-1)+1]=0; rampright[4*(i*ts+ts-1)+2]=0; rampright[4*(i*ts+ts-1)+3]=0;

            preinttab[4*((ts-1)*ts+i)]=0; preinttab[4*((ts-1)*ts+i)+1]=0; preinttab[4*((ts-1)*ts+i)+2]=0; preinttab[4*((ts-1)*ts+i)+3]=0;
            rampleft[4*((ts-1)*ts+i)]=0; rampleft[4*((ts-1)*ts+i)+1]=0; rampleft[4*((ts-1)*ts+i)+2]=0; rampleft[4*((ts-1)*ts+i)+3]=0;
            rampright[4*((ts-1)*ts+i)]=0; rampright[4*((ts-1)*ts+i)+1]=0; rampright[4*((ts-1)*ts+i)+2]=0; rampright[4*((ts-1)*ts+i)+3]=0;
        }
    }
};

// DATA for EXAMPLE RENDERING

int W=1000;
int H=800;

Raw3DData obj;
unsigned char *tf;

float mm[16];
mat3<float> mat;
vec3<float> vi(1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0),
vj,
vk(-1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0));

int mov=0;
int posx, posy;
float alpha=0.0, beta=0.0;

unsigned int vshader;
const int HV_VOLDATA_NVS = 21;
const char *hvVolDataVS[HV_VOLDATA_NVS] = {
    "#version 140\n",
    "\n",
    "in vec4 vpos;\n",
    "uniform mat4 trans;\n",
    "uniform mat4 rev;\n",
    "uniform mat4 persp;\n",
    "uniform vec3 eye;\n",
    "uniform vec3 light;\n",
    "\n",
    "smooth out vec3 textcoord;\n",
    "noperspective out vec3 dir;\n",
    "noperspective out vec3 lidir;\n",
    //"smooth out vec3 view;\n",
    "\n",
    "void main(void) {\n",
    " vec4 v = trans * vpos;\n",
    " gl_Position = persp * v;\n",
    " vec3 view = v.xyz/v.w-eye;\n",
    " dir = normalize(mat3(rev) * view);\n",
    " lidir = normalize(mat3(rev) * light);\n",
    " textcoord = vec3((vpos.x+1.0)/2.0,(vpos.y+1.0)/2.0,(vpos.z+1.0)/2.0);\n",
    "}\n"
};

unsigned int fshader;
// USING NO PRE-INTEGRATION
const int HV_VOLDATA_NFS = 31;
const char *hvVolDataFS[HV_VOLDATA_NFS] = {
    "#version 140\n",
    "\n",
    "smooth in vec3 textcoord;\n",
    "noperspective in vec3 dir;\n",
    "out vec3 gl_FragColor;\n",
    "uniform sampler1D tf;\n",
    "uniform sampler3D data;\n",
    "noperspective in vec3 lidir;\n",
    "uniform vec4 lipara;\n",
    "uniform float slice, sampling;\n",
    "void main(void) {\n",
    " vec3 ndir = normalize(dir);\n",
    " vec3 ldir = normalize(lidir);\n",
    " vec3 dstep = ndir*(sampling/128.0);\n",
    " vec3 pos = textcoord+dstep;\n",
    " vec4 col=vec4(0.0,0.0,0.0,1.0);\n",
    " while(pos.x<1.0 && pos.x>0.0 && pos.y<1.0 && pos.y>0.0 && pos.z<1.0 && pos.z>0.0  && col.w>0.01)\n",
    "  {\n",
    "   vec4 val = vec4(texture(data, pos));\n",
    "   vec4 cc = vec4(texture(tf, val.w));\n",
    "   vec3 grad = normalize(vec3(val.x*2.0-1.0, val.y*2.0-1.0, val.z*2.0-1.0));\n",
    "	float alpha = col.w*cc.w*sampling;\n",
    "	vec3 half; if (dot(grad,ldir)<0.0 || dot(grad,ndir)>=0.0) half=vec3(0.0);\n",
    "	else half = normalize(ldir-ndir);\n",
    "	float dd = clamp(dot(ldir,grad),0.0,1.0); \n",
    "	float ds = pow(clamp(dot(grad,half),0.0,1.0),lipara.w);\n",
    "	col = vec4(col.xyz+(cc.xyz*(lipara.x+lipara.y*dd)+vec3(ds*lipara.z))*alpha,col.w*(1.0-cc.w*sampling));\n",
    "   pos = pos+dstep;\n",
    "  }\n",
    " gl_FragColor = col.xyz+vec3(col.w);\n",
    "}"
};

// USING CLASSIC  PRE-INTEGRATION
const int HV_VOLDATA_NFS2 = 35;
const char *hvVolDataFS2[HV_VOLDATA_NFS2] = {
    "#version 140\n",
    "\n",
    "smooth in vec3 textcoord;\n",
    "noperspective in vec3 dir;\n",
    "out vec3 gl_FragColor;\n",
    "uniform sampler2D pretab;\n",
    "uniform sampler3D data;\n",
    "noperspective in vec3 lidir;\n",
    "uniform vec4 lipara;\n",
    "\n",
    "void main(void) {\n",
    " vec3 ndir = normalize(dir);\n",
    " vec3 ldir = normalize(lidir); \n",
    " vec3 step = ndir*(0.5/128.0);\n",
    " vec3 pos = textcoord+step;\n",
    " vec4 col=vec4(0.0,0.0,0.0,1.0);\n",
    " vec4 oldval = vec4(texture(data, pos));\n",
    " vec3 oldgrad = normalize(vec3(oldval.x*2.0-1.0, oldval.y*2.0-1.0, oldval.z*2.0-1.0));\n",
    " pos = pos+step*0.5;\n",
    " while(pos.x<1.0 && pos.x>0.0 && pos.y<1.0 && pos.y>0.0 && pos.z<1.0 && pos.z>0.0 && col.w>0.01)\n",
    "  {\n",
    "   vec4 val = vec4(texture(data, pos));\n",
    "   vec4 cc = vec4(texture(pretab, vec2(oldval.w, val.w)));\n",
    "	float alpha = col.w;\n",
    "   vec3 grad = normalize(vec3(val.x*2.0-1.0, val.y*2.0-1.0, val.z*2.0-1.0));\n",
    "	vec3 half; if (dot(grad,ldir)<0.0 || dot(grad,ndir)>=0.0) half=vec3(0.0);\n",
    "	else half = normalize(ldir-ndir);\n",
    "	float dd = clamp(dot(ldir,grad),0.0,1.0); \n",
    "	float ds = clamp((dot(grad,half)-1.0)*lipara.w+1.0,0.0,1.0);\n",
    "	col = vec4(col.xyz+(cc.xyz*(lipara.x+lipara.y*dd)+vec3(ds*lipara.z))*alpha,alpha*(1.0-cc.w));\n",
    "   pos = pos+step;\n",
    "	oldval = val; oldgrad=grad; \n",
    "  }\n",
    " gl_FragColor = col.xyz+vec3(col.w);\n",
    "}"
};

unsigned int id[4], textid[4];
unsigned int gpuprog;
int itrans, iposition, ipersp, ilight, ieye, ipara, islice, isampling;
int isampler, itf;

float vert[3*8] = { -1.0f, -1.0f, -1.0f,   -1.0f, -1.0f, 1.0f,   -1.0f, 1.0f, -1.0f,   -1.0f, 1.0f, 1.0f,
                    1.0f, -1.0f, -1.0f,   1.0f, -1.0f, 1.0f,   1.0f, 1.0f, -1.0f,   1.0f, 1.0f, 1.0f };
//float vert[3*4] = { 0.0f,-1.0f,-1.0f,0.0f,-1.0f,1.0f, 0.0f,1.0f,1.0f, 0.0f,1.0f,-1.0f };
unsigned int faces[36] = { 0,1,2, 2,1,3 , 3,1,5, 3,5,7,   0,2,4, 4,2,6, 1,0,4, 1,4,5, 2,3,7, 2,7,6, 6,7,5, 5,4,6 };

void createlist()
{
    int i;
    int compiled=0;
    char buffer[1024];

    // create buffers for vertices, normals and element indices
    glGenBuffers(4,id);
    glBindBuffer(GL_ARRAY_BUFFER, id[0]);
    glBufferData(GL_ARRAY_BUFFER, 3*8*sizeof (float), vert, GL_STATIC_DRAW);
    CHECK_GL_ERROR();

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*12*sizeof (unsigned int), faces, GL_STATIC_DRAW);
    CHECK_GL_ERROR();

    // vertex shader
    int vshaderlength[500];
    vshader = glCreateShader(GL_VERTEX_SHADER);
    for (i=0; i<HV_VOLDATA_NVS; i++) vshaderlength[i]=strlen(hvVolDataVS[i]);
    glShaderSource(vshader, HV_VOLDATA_NVS, (const GLchar **)hvVolDataVS, vshaderlength);
    glCompileShader(vshader);
    CHECK_GL_ERROR();
    glGetShaderiv(vshader, GL_COMPILE_STATUS, &compiled);
    CHECK_GL_ERROR();
    if (compiled) { printf("Vertex Shader compilation successful!\n"); }
    else
    {
        printf("Shader compile error:\n");
        glGetShaderInfoLog(vshader, 1024, &i, buffer);
        printf("%s\n", buffer);
    }
    // fragment shader
    int fshaderlength[500];
    fshader = glCreateShader(GL_FRAGMENT_SHADER);
#ifdef PREINT
    for (i=0; i<HV_VOLDATA_NFS2; i++) fshaderlength[i]=strlen(hvVolDataFS2[i]);
    glShaderSource(fshader, HV_VOLDATA_NFS2, (const GLchar **)hvVolDataFS2, fshaderlength);
#else
    for (i=0; i<HV_VOLDATA_NFS; i++) fshaderlength[i]=strlen(hvVolDataFS[i]);
    glShaderSource(fshader, HV_VOLDATA_NFS, (const GLchar **)hvVolDataFS, fshaderlength);
#endif
    glCompileShader(fshader);
    CHECK_GL_ERROR();
    glGetShaderiv(fshader, GL_COMPILE_STATUS, &compiled);
    CHECK_GL_ERROR();
    if (compiled) { printf("Fragment Shader compilation successful!\n"); }
    else
    {
        printf("Shader compile error:\n");
        glGetShaderInfoLog(fshader, 1024, &i, buffer);
        printf("%s\n", buffer);
    }

    gpuprog = glCreateProgram();
    glAttachShader(gpuprog, vshader);
    glAttachShader(gpuprog, fshader);
    glLinkProgram(gpuprog);
    CHECK_GL_ERROR();
    glGetProgramiv(gpuprog, GL_LINK_STATUS, &compiled);
    CHECK_GL_ERROR();
    if (compiled) { printf("Link successful!\n"); }
    else
    {
        printf("Link error:\n");
        glGetProgramInfoLog(gpuprog, 1024, &i, buffer);
        printf("%s\n", buffer);
    }

    glUseProgram(gpuprog);
    CHECK_GL_ERROR();

    int nattr; glGetProgramiv(gpuprog, GL_ACTIVE_ATTRIBUTES, &nattr);
    printf("activated attributes : %d\n", nattr);
    int nunif; glGetProgramiv(gpuprog, GL_ACTIVE_UNIFORMS, &nunif);
    printf("activated uniforms : %d\n", nunif);

    int ivpos = glGetAttribLocation(gpuprog, "vpos");
    printf("attribute vpos = %d\n", ivpos);
    CHECK_GL_ERROR();

    itrans = glGetUniformLocation(gpuprog, "trans");
    CHECK_GL_ERROR();
    printf("uniform trans = %d\n", itrans);
    iposition = glGetUniformLocation(gpuprog, "rev");
    CHECK_GL_ERROR();
    printf("uniform rev = %d\n", iposition);
    ipersp = glGetUniformLocation(gpuprog, "persp");
    CHECK_GL_ERROR();
    printf("uniform persp = %d\n", ipersp);
    ieye = glGetUniformLocation(gpuprog, "eye");
    CHECK_GL_ERROR();
    printf("uniform eyepos = %d\n", ieye);
    ilight = glGetUniformLocation(gpuprog, "light");
    CHECK_GL_ERROR();
    printf("uniform light = %d\n", ilight);
    ipara = glGetUniformLocation(gpuprog, "lipara");
    CHECK_GL_ERROR();
    printf("uniform light = %d\n", ilight);
    isampler = glGetUniformLocation(gpuprog, "data");
    CHECK_GL_ERROR();
    printf("sampler3D data = %d\n", isampler);

#ifdef PREINT
    itf = glGetUniformLocation(gpuprog, "pretab");
    CHECK_GL_ERROR();
    printf("sampler1D pretab = %d\n", itf);
#else
    islice = glGetUniformLocation(gpuprog, "slice");
    CHECK_GL_ERROR();
    printf("uniform slice = %d\n", islice);
    isampling = glGetUniformLocation(gpuprog, "sampling");
    CHECK_GL_ERROR();
    printf("uniform sampling = %d\n", isampling);
    itf = glGetUniformLocation(gpuprog, "tf");
    CHECK_GL_ERROR();
    printf("sampler1D tf = %d\n", itf);
#endif

    glBindBuffer(GL_ARRAY_BUFFER, id[0]);
    glVertexAttribPointer(ivpos, 3, GL_FLOAT, false, 0, NULL);
    CHECK_GL_ERROR();
    glEnableVertexAttribArray(ivpos);
    CHECK_GL_ERROR();
}

void inittextures()
{ 
    glGenTextures(4, textid);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, textid[0]);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA8, obj.sx, obj.sy, obj.sz, 0, GL_RGBA, GL_UNSIGNED_BYTE, obj.grad);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

#ifdef PREINT
    glActiveTexture(GL_TEXTURE0+1);
    glBindTexture(GL_TEXTURE_2D, textid[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16, obj.tabsize, obj.tabsize, 0, GL_RGBA, GL_UNSIGNED_SHORT, obj.preinttab);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
#else
    glActiveTexture(GL_TEXTURE0+1);
    glBindTexture(GL_TEXTURE_1D, textid[1]);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA8, obj.tabsize, 0, GL_RGBA, GL_UNSIGNED_BYTE, obj.transfer);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
#endif
}

// OPENGL DISPLAY
void display()
{
    transform<float> tr;

    // effacer l'\E9cran
    glClearColor(1.0,1.0,1.0,0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );

    // initialiser le viewport
    glViewport(0, 0, W, H);

    // activer le Z-buffer
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glFrontFace(GL_CCW);
    glEnable(GL_CULL_FACE);


    float li[3] = { 1.0, 1.0, 1.0 };
    glUniform3fv(ilight, 1, li);
    CHECK_GL_ERROR();
    float liparadata[4];
    // define Ka, Kd,Ks, Phong exp
    liparadata[0]=0.3; liparadata[1]=0.7; liparadata[2]=0.1; liparadata[3]=60.0;
    glUniform4fv(ipara, 1, liparadata);
    CHECK_GL_ERROR();

#ifdef PREINT
#else
    float slice_val=0.0;
    glUniform1fv(islice, 1, &slice_val);
    float sampling_val= 0.5 ; //voxelstep;
    glUniform1fv(isampling, 1, &sampling_val);
#endif
    // initialiser la transformation perspective
    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glFrustum(-0.42, 0.42, -0.42*H/W, 0.42*H/W, 1.0, 50.0);


    // initialiser la transformation de l'objet
    mat = mat3<float>(vi,vj,vk);
    mat.inverse();
    mm[0]=mat.I().X(); mm[1]=mat.I().Y(); mm[2]=mat.I().Z(); mm[3]=0.0;
    mm[4]=mat.J().X(); mm[5]=mat.J().Y(); mm[6]=mat.J().Z(); mm[7]=0.0;
    mm[8]=mat.K().X(); mm[9]=mat.K().Y(); mm[10]=mat.K().Z(); mm[11]=0.0;
    mm[12]=0.0; mm[13]=0.0; mm[14]=0.0; mm[15]=1.0;

    float eye[3]={ -12.0, -12.0, 12.0 };
    glUniform3fv(ieye, 1, eye);
    CHECK_GL_ERROR();

    tr.setIdentity();
    tr.frustum(-0.1, 0.1, -0.1*H/W, 0.1*H/W, 0.08, 50.0);
    tr.frame(vi,vj,vk);
    tr.translation(vec3<float>(-eye[0], -eye[1], -eye[2]));
    tr.getMatrix(mm);
    glUniformMatrix4fv(ipersp, 1, false, mm);
    CHECK_GL_ERROR();

    tr.setIdentity();
    tr.rotation(vec3<float>(0.0,0.0,1.0), alpha/180.0*PI);
    tr.rotation(vec3<float>(1.0,0.0,0.0), beta/180.0*PI);
    tr.scale(vec3<float>(8.0,8.0,8.0));
    tr.getMatrix(mm);
    glUniformMatrix4fv(itrans, 1, false, mm);
    CHECK_GL_ERROR();
    tr.setIdentity();
    tr.scale(vec3<float>(1.0/8.0,1.0/8.0,1.0/8.0));
    tr.rotation(vec3<float>(1.0,0.0,0.0), -beta/180.0*PI);
    tr.rotation(vec3<float>(0.0,0.0,1.0), -alpha/180.0*PI);
    tr.getMatrix(mm);
    glUniformMatrix4fv(iposition, 1, false, mm);
    CHECK_GL_ERROR();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, textid[0]);
    CHECK_GL_ERROR();
    glUniform1i(isampler,0);
    CHECK_GL_ERROR();
    glEnable(GL_TEXTURE_3D);
    CHECK_GL_ERROR();
    glActiveTexture(GL_TEXTURE0+1);
#ifdef PREINT
    glBindTexture(GL_TEXTURE_2D, textid[1]);
    CHECK_GL_ERROR();
    glUniform1i(itf,1);
    CHECK_GL_ERROR();
    glEnable(GL_TEXTURE_2D);
    CHECK_GL_ERROR();
#else
    glBindTexture(GL_TEXTURE_1D, textid[1]);
    CHECK_GL_ERROR();
    glUniform1i(itf,1);
    CHECK_GL_ERROR();
    glEnable(GL_TEXTURE_1D);
    CHECK_GL_ERROR();
#endif

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id[1]);
    glDrawElements(GL_TRIANGLES, 12*3, GL_UNSIGNED_INT, NULL);
    CHECK_GL_ERROR();

    glutSwapBuffers();
}


// MANAGE EVENTS
void clic(int button, int state, int x, int y)
{
    if (button==0 && state==0) { mov=1; posx=x; posy=y; }
    else mov=0;
}

void bouge(int x, int y)
{
    if (mov==1)
    {
        //printf("move %d,%d\n", x,y);
        alpha += (float)(x-posx)*0.5;
        posx=x;
        beta += (float)(y-posy)*0.25;
        posy=y;
        display();
    }
}

void key(unsigned char c, int x, int y)
{
    printf("char=%c (%d)\n",c, (int)c);
}



// START PROGRAMM AND INITIALIZE WINDOW
int main(int argc, char **argv)
{

    vj.cross(vi,vk);
    vj.reverse();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutInitWindowPosition(0,0) ;
    glutInitWindowSize(W, H);
    int win = glutCreateWindow("ma fenetre");

    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }

    glutDisplayFunc(display);
    glutMouseFunc(clic);
    glutMotionFunc(bouge);
    //  glutKeyboardFunc(key);

    // fonction de transfert
    tf = new unsigned char [128*4];
    int i;
    for (i=0; i<127; i++) { tf[4*i]=0; tf[4*i+1]=0; tf[4*i+2]=0; tf[4*i+3]=0; }
    for (i=10; i<30; i++) { tf[4*i]=static_cast<unsigned char>(0.99*255.0); tf[4*i+1]=(unsigned char)(0.9*255.0); tf[4*i+2]=(unsigned char)(0.6*255.0); tf[4*i+3]=(unsigned char)(0.025*255.0); }
    for (i=30; i<40; i++) { tf[4*i]=(unsigned char)(0.99*255.0); tf[4*i+1]=(unsigned char)(0.29*255.0); tf[4*i+2]=(unsigned char)(0.1*255.0); tf[4*i+3]=(unsigned char)(0.2*255.0); }
    for (i=40; i<50; i++) { tf[4*i]=(unsigned char)(0.5*255.0); tf[4*i+1]=(unsigned char)(0.7*255.0); tf[4*i+2]=(unsigned char)(0.9*255.0); tf[4*i+3]=(unsigned char)(0.85*255.0); }
    for (i=70; i<90; i++) { tf[4*i]=(unsigned char)(0.99*255.0); tf[4*i+1]=(unsigned char)(0.8*255.0); tf[4*i+2]=(unsigned char)(0.7*255.0); tf[4*i+3]=(unsigned char)(0.85*255.0); }
    for (i=100; i<126; i++) { tf[4*i]=(unsigned char)(0.99*255.0); tf[4*i+1]=(unsigned char)(0.9*255.0); tf[4*i+2]=(unsigned char)(0.6*255.0); tf[4*i+3]=(unsigned char)(0.85*255.0); }

    // charger les donnees
    printf("loading raw data...\n");
    obj.importRAW("data/CTKnee.raw",379,229,305, 1.0/255.0);
    //obj.computeGrad();
    obj.importGrad("data/CTKnee_grad.raw");
    obj.computePreIntegration(tf, 128, 0.5, 20);

    // init GPU
    createlist();
    inittextures();

    glutMainLoop();
    return 0;
}

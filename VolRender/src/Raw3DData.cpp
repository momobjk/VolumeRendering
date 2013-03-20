#include "Raw3DData.h"

Raw3DData::Raw3DData()
{
    sx=0; sy=0; sz=0;
    data=0; grad=0; scale = 1.0/255.0;
}

unsigned char Raw3DData::get(int ind[3]) const
{
    return data[ind[0]+ind[1]*sx+ind[2]*sx*sy];
}

unsigned char Raw3DData::get(int i, int j, int k) const
{
    return data[i+j*sx+k*sx*sy];
}

double Raw3DData::rawdensity(int i, int j, int k, double ri, double rj, double rk) const
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

vec3<double> Raw3DData::rawgradient(int i, int j, int k, double ri, double rj, double rk) const
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
bool Raw3DData::importRAW(const char *name, int xx, int yy, int zz, double ss)
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

bool Raw3DData::importGrad(const char *name)
{

    FILE *fd = fopen(name, "rb");
    if (fd==0) { printf("cannot open file!\n"); exit(1); }
    grad = new unsigned char [4*sx*sy*sz];
    fread(grad,sx*sy*sz*4, sizeof(unsigned char), fd);
    fclose(fd);
    return true;
}

void Raw3DData::computeGrad()
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

void Raw3DData::computePreIntegration(unsigned char *tf, int ts, double vs, int nstep)
    {
        int i,j,k;
        double scal=65535.0;

        tabSize=ts;
        voxelStep=vs;
        transfer = tf;

        preIntTab=(unsigned short *)malloc(ts*ts*4*sizeof (unsigned short));
        rampLeft=(unsigned short *)malloc(ts*ts*4*sizeof (unsigned short));
        rampRight=(unsigned short *)malloc(ts*ts*4*sizeof (unsigned short));
        for (i=0; i<tabSize*tabSize; i++) { preIntTab[i]=0; rampLeft[i]=0; rampRight[0]=0; }

        for (i=0; i<tabSize; i++)
        {
            double alpha=0.0, colt=1.0, alphamid=0.0;
            double vr = 0.0, vg = 0.0, vb = 0.0, alpha0r=0.0;
            double vrl = 0.0, vgl = 0.0, vbl = 0.0, alphal=0.0;
            double vrr = 0.0, vgr = 0.0, vbr = 0.0, alphar=0.0;
            double vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
            double vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;

            double colr, colg,colb, de;
            //printf("step %d/%d\n", i,pi.sizeX());
            alpha=0.0; colt=1.0; alphamid=0.0; vr = 0.0; vg = 0.0; vb = 0.0; alpha0r=0.0;
            vrl = 0.0; vgl = 0.0; vbl = 0.0; alphal=0.0;
            vrr = 0.0; vgr = 0.0; vbr = 0.0; alphar=0.0;
            vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
            vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;
            colr = (double)tf[4*i+0]/255.0;
            colg = (double)tf[4*i+1]/255.0;
            colb = (double)tf[4*i+2]/255.0;
            de = (double)tf[4*i+3]/255.0;

            for (k=0; k<nstep; k++)
            {
                double tt = (double)k/(double)(nstep-1);
                alpha=colt*de/(double)nstep*voxelStep;
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
                colt *= (1.0-de/(double)nstep*voxelStep);
            }
            preIntTab[4*(i*ts+i)]=(unsigned short)(scal*vr); preIntTab[4*(i*ts+i)+1]=(unsigned short)(scal*vg); preIntTab[4*(i*ts+i)+2]=(unsigned short)(scal*vb); preIntTab[4*(i*ts+i)+3]=(unsigned short)(scal*(1.0-colt));
            rampLeft[4*(i*ts+i)]=(unsigned short)(scal*vrl); rampLeft[4*(i*ts+i)+1]=(unsigned short)(scal*vgl); rampLeft[4*(i*ts+i)+2]=(unsigned short)(scal*vbl); rampLeft[4*(i*ts+i)+3]=(unsigned short)(scal*alphal);
            rampRight[4*(i*ts+i)]=(unsigned short)(scal*vrr); rampRight[4*(i*ts+i)+1]=(unsigned short)(scal*vgr); rampRight[4*(i*ts+i)+2]=(unsigned short)(scal*vbr); rampRight[4*(i*ts+i)+3]=(unsigned short)(scal*alphar);
            printf("i=%d  Ka=%g,%g,%g   Kd=%g,%g,%g\n",i,vr,vg,vb,vrr,vgr,vbr);
            for (j=i+1; j<ts; j++)
            {
                double d = 1.0/(double)(j-i);
                //double kk=0.0;
                alpha=0.0; colt=1.0; alphamid=0.0;  vr = 0.0; vg = 0.0; vb = 0.0; alpha0r=0.0;
                vrl = 0.0; vgl = 0.0; vbl = 0.0; alphal=0.0;
                vrr = 0.0; vgr = 0.0; vbr = 0.0; alphar=0.0;
                vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
                vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;

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

                    alpha = colt*de/(double)nstep*d*voxelStep;
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
                    colt *= (1.0-de/(double)nstep*d*voxelStep);

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
                preIntTab[4*(j*ts+i)]=(unsigned short)(scal*vr); preIntTab[4*(j*ts+i)+1]=(unsigned short)(scal*vg); preIntTab[4*(j*ts+i)+2]=(unsigned short)(scal*vb); preIntTab[4*(j*ts+i)+3]=(unsigned short)(scal*(1.0-colt));
                rampLeft[4*(j*ts+i)]=(unsigned short)(scal*vrl); rampLeft[4*(j*ts+i)+1]=(unsigned short)(scal*vgl); rampLeft[4*(j*ts+i)+2]=(unsigned short)(scal*vbl); rampLeft[4*(j*ts+i)+3]=(unsigned short)(scal*alphal);
                rampRight[4*(j*ts+i)]=(unsigned short)(scal*vrr); rampRight[4*(j*ts+i)+1]=(unsigned short)(scal*vgr); rampRight[4*(j*ts+i)+2]=(unsigned short)(scal*vbr); rampRight[4*(j*ts+i)+3]=(unsigned short)(scal*alphar);
            }
            for (j=0; j<i; j++)
            {
                double d = 1.0/(double)(i-j);
                //double kk=0.0;
                colt=1.0; alpha=0.0; alphamid=0.0; vr = 0.0; vg = 0.0; vb = 0.0; alpha0r=0.0;
                vrl = 0.0; vgl = 0.0; vbl = 0.0; alphal=0.0;
                vrr = 0.0; vgr = 0.0; vbr = 0.0; alphar=0.0;
                vr2r = 0.0, vg2r = 0.0, vb2r = 0.0, alpha2r=0.0;
                vr3r = 0.0, vg3r = 0.0, vb3r = 0.0, alpha3r=0.0;

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

                    alpha = colt*de/(double)nstep*d*voxelStep;
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
                    colt *= (1.0-de/(double)nstep*d*voxelStep);
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
                preIntTab[4*(j*ts+i)]=(unsigned short)(scal*vr); preIntTab[4*(j*ts+i)+1]=(unsigned short)(scal*vg); preIntTab[4*(j*ts+i)+2]=(unsigned short)(scal*vb); preIntTab[4*(j*ts+i)+3]=(unsigned short)(scal*(1.0-colt));
                rampLeft[4*(j*ts+i)]=(unsigned short)(scal*vrl); rampLeft[4*(j*ts+i)+1]=(unsigned short)(scal*vgl); rampLeft[4*(j*ts+i)+2]=(unsigned short)(scal*vbl); rampLeft[4*(j*ts+i)+3]=(unsigned short)(scal*alphal);
                rampRight[4*(j*ts+i)]=(unsigned short)(scal*vrr); rampRight[4*(j*ts+i)+1]=(unsigned short)(scal*vgr); rampRight[4*(j*ts+i)+2]=(unsigned short)(scal*vbr); rampRight[4*(j*ts+i)+3]=(unsigned short)(scal*alphar);
            }
        }
        for (i=0; i<ts; i++)
        {
            preIntTab[4*(i*ts+0)]=0; preIntTab[4*(i*ts+0)+1]=0; preIntTab[4*(i*ts+0)+2]=0; preIntTab[4*(i*ts+0)+3]=0;
            rampLeft[4*(i*ts+0)]=0; rampLeft[4*(i*ts+0)+1]=0; rampLeft[4*(i*ts+0)+2]=0; rampLeft[4*(i*ts+0)+3]=0;
            rampRight[4*(i*ts+0)]=0; rampRight[4*(i*ts+0)+1]=0; rampRight[4*(i*ts+0)+2]=0; rampRight[4*(i*ts+0)+3]=0;

            preIntTab[4*(0*ts+i)]=0; preIntTab[4*(0*ts+i)+1]=0; preIntTab[4*(0*ts+i)+2]=0; preIntTab[4*(0*ts+i)+3]=0;
            rampLeft[4*(0*ts+i)]=0; rampLeft[4*(0*ts+i)+1]=0; rampLeft[4*(0*ts+i)+2]=0; rampLeft[4*(0*ts+i)+3]=0;
            rampRight[4*(0*ts+i)]=0; rampRight[4*(0*ts+i)+1]=0; rampRight[4*(0*ts+i)+2]=0; rampRight[4*(0*ts+i)+3]=0;

            preIntTab[4*(i*ts+ts-1)]=0; preIntTab[4*(i*ts+ts-1)+1]=0; preIntTab[4*(i*ts+ts-1)+2]=0; preIntTab[4*(i*ts+ts-1)+3]=0;
            rampLeft[4*(i*ts+ts-1)]=0; rampLeft[4*(i*ts+ts-1)+1]=0; rampLeft[4*(i*ts+ts-1)+2]=0; rampLeft[4*(i*ts+ts-1)+3]=0;
            rampRight[4*(i*ts+ts-1)]=0; rampRight[4*(i*ts+ts-1)+1]=0; rampRight[4*(i*ts+ts-1)+2]=0; rampRight[4*(i*ts+ts-1)+3]=0;

            preIntTab[4*((ts-1)*ts+i)]=0; preIntTab[4*((ts-1)*ts+i)+1]=0; preIntTab[4*((ts-1)*ts+i)+2]=0; preIntTab[4*((ts-1)*ts+i)+3]=0;
            rampLeft[4*((ts-1)*ts+i)]=0; rampLeft[4*((ts-1)*ts+i)+1]=0; rampLeft[4*((ts-1)*ts+i)+2]=0; rampLeft[4*((ts-1)*ts+i)+3]=0;
            rampRight[4*((ts-1)*ts+i)]=0; rampRight[4*((ts-1)*ts+i)+1]=0; rampRight[4*((ts-1)*ts+i)+2]=0; rampRight[4*((ts-1)*ts+i)+3]=0;
        }
    }

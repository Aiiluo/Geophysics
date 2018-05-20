/* 2-D two-components wavefield modeling using pseudo-pure mode P-wave equation in VTI media.
    Jiubing Cheng and Wei Kang 

                        modify: 2018.02
*/
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#ifndef pi
#define pi 3.141592653
#endif

#ifdef _OPENMP
#include<omp.h>
#endif

float WltRicker(float t, float favg, float Amp) 
{
        float tdelay = 1.0/favg;
        float x=pow(pi*favg*(t-tdelay),2);
        return -Amp*exp(-x)*(1-2*x);
}

void zero1float(float *p, int n1)
/*< free a 1-d array of float >*/
{
     int i;
	 #ifdef _OPENMP
	 #pragma omp parallel for private(i) 
	 #endif
     for(i=0;i<n1;i++) p[i]=0.0;
}

void zero2float(float **p, int n1, int n2)
/*< free a 2-d array of float >*/
{
     int i, j;
	 #ifdef _OPENMP
	 #pragma omp parallel for private(i,j) 
	 #endif
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=0.0;
}

/*************************************************************************
 * Calculate phase velocity, group angle amd group velocity for
 * qP, qSV-waves in 2D VTI media
 *************************************************************************/
#ifndef pi
#define pi 3.141592653
#endif
float VPphase (float vp0, float vs0, float epsilon, float delta, float ap)
/*< VPphase : calculate phase velocity for given  
 * phase angle of qP-wave in 2D VTI media >*/
/* Tsvanki */
/* An acoustic wave equation for pure P wave int 2D-TTI media, 2011SEG
    Paul L. Stoffa, University of Texas at Austin */
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;
    float vpp;

    vpp = vp0*sqrt(1+epsilon*sin(ap)*sin(ap)
        - f/2.0
        + f/2.0*sqrt(pow(1+2*epsilon*sin(ap)*sin(ap)/f,2)
                    -2*(epsilon-delta)*sin(2*ap)*sin(2*ap)/f));

    return vpp;
}

float VSphase (float vp0, float vs0, float epsilon, float delta, float ap)
/*< VSphase : calculate phase velocity for given 
 * phase angle of qSV-wave in 2D VTI media >*/
/* ap:Symmetry axis phase angle */
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;
    float vsp;

    vsp = vp0*sqrt(1+epsilon*sin(ap)*sin(ap)
        - f/2.0
        - f/2.0*sqrt(pow(1+2*epsilon*sin(ap)*sin(ap)/f,2)
                    -2*(epsilon-delta)*sin(2*ap)*sin(2*ap)/f));

    return vsp;
}

float VPgroup (float vp0, float vs0, float epsilon, float delta, float ap)
/*< VPgroup : calculate group velocity from 
 * phase velocity of qP-wave in 2D VTI media >*/
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;

    float vpp = VPphase (vp0, vs0, epsilon, delta, ap);
        
    float A = epsilon*sin(2*ap)
            + f*f/4.0/sqrt(4*epsilon*epsilon*pow(sin(ap),4)
                + f*4.0*sin(ap)*sin(ap)*(2*delta*cos(ap)*cos(ap)-epsilon*cos(2*ap))+f*f)
            * (4*epsilon*sin(2*ap)/f+8*epsilon*epsilon*sin(ap)*sin(ap)*sin(2*ap)/(f*f)
                - 8*(epsilon-delta)*sin(2*ap)*cos(2*ap)/f);

    float dvda = A*vp0*vp0/vpp/2.0;

    float vpg = vpp*sqrt(1+pow(dvda/vpp,2));
    return vpg;
}

float VSgroup (float vp0, float vs0, float epsilon, float delta, float ap)
/*< VSgroup : calculate group velocity from 
 * phase velocity of qSV-wave in 2D VTI media >*/
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;

    float vsp = VSphase (vp0, vs0, epsilon, delta, ap);
        
    float A = epsilon*sin(2*ap)
                - f*f/4.0/sqrt(4*epsilon*epsilon*pow(sin(ap),4)
                    + f*4.0*sin(ap)*sin(ap)*(2*delta*cos(ap)*cos(ap)-epsilon*cos(2*ap))+f*f)
                * (4*epsilon*sin(2*ap)/f+8*epsilon*epsilon*sin(ap)*sin(ap)*sin(2*ap)/(f*f)
                    - 8*(epsilon-delta)*sin(2*ap)*cos(2*ap)/f);

    float dvda = A*vp0*vp0/vsp/2.0;

    float vsg=vsp*sqrt(1+pow(dvda/vsp,2));

    return vsg;
}

float AnglePgroup (float vp0, float vs0, float epsilon, float delta, float ap)
/*< AnglePgroup : calculate group angle from phase velocity 
 * and angle of qP-wave in 2D VTI media >*/
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;

    float vpp = VPphase (vp0, vs0, epsilon, delta, ap);

    float A = epsilon*sin(2*ap)
                + f*f/4.0/sqrt(4*epsilon*epsilon*pow(sin(ap),4)
                    + f*4.0*sin(ap)*sin(ap)*(2*delta*cos(ap)*cos(ap)-epsilon*cos(2*ap))+f*f)
                * (4*epsilon*sin(2*ap)/f+8*epsilon*epsilon*sin(ap)*sin(ap)*sin(2*ap)/(f*f)
                    - 8*(epsilon-delta)*sin(2*ap)*cos(2*ap)/f);

    float dvda = A*vp0*vp0/vpp/2.0;

    float B = 1-tan(ap)*dvda/vpp;
    float D = tan(ap)+dvda/vpp;

    float C, apg;

    if(fabs(B)<1e-30)
        B=1e-30;

    C=D/B;
    if(fabs(C)>1e30)
        C=1e30;

    apg=atan(C);
    if(apg<0)
        apg+=pi;
	
    return apg;
}

float AngleSgroup (float vp0, float vs0, float epsilon, float delta, float ap)
/*< AngleSgroup : calculate group angle from phase velocity 
 * and angle of qSV-wave in 2D VTI media >*/
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;

    float vsp = VSphase (vp0, vs0, epsilon, delta, ap);

    float A = epsilon*sin(2*ap)
                - f*f/4.0/sqrt(4*epsilon*epsilon*pow(sin(ap),4)
                    + f*4.0*sin(ap)*sin(ap)*(2*delta*cos(ap)*cos(ap)-epsilon*cos(2*ap))+f*f)
                * (4*epsilon*sin(2*ap)/f+8*epsilon*epsilon*sin(ap)*sin(ap)*sin(2*ap)/(f*f)
                    - 8*(epsilon-delta)*sin(2*ap)*cos(2*ap)/f);

    float dvda = A*vp0*vp0/vsp/2.0;

    float B=1-tan(ap)*dvda/vsp;
    float D=tan(ap)+dvda/vsp;

    float C, asg;

    if(fabs(B)<1e-20)
        B=1e-20;

    C=D/B;
    if(fabs(C)>1e5)
        C=1e5;

    asg=atan(C);
    if(asg<0)
        asg+=pi;
	
    return asg;
}

void VPAnglePgroup (float vp0, float vs0, float epsilon, float delta, float ap,
                   float *vpg, float *apg)
/*< VPAnglePgroup : calculate group velocity and angle from phase velocity for
 * given phase angle of qP-wave in 2D VTI media >*/
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;

    float vpp = VPphase (vp0, vs0, epsilon, delta, ap);
        
    float A = epsilon*sin(2*ap)
                + f*f/4.0/sqrt(4*epsilon*epsilon*pow(sin(ap),4)
                    + f*4.0*sin(ap)*sin(ap)*(2*delta*cos(ap)*cos(ap)-epsilon*cos(2*ap))+f*f)
                * (4*epsilon*sin(2*ap)/f+8*epsilon*epsilon*sin(ap)*sin(ap)*sin(2*ap)/(f*f)
                    - 8*(epsilon-delta)*sin(2*ap)*cos(2*ap)/f);

    float dvda = A*vp0*vp0/vpp/2.0;

    float B, D, C;

    *vpg = vpp*sqrt(1+pow(dvda/vpp,2));

    B = 1-tan(ap)*dvda/vpp;
    D = tan(ap)+dvda/vpp;

    if(fabs(B)<1e-30)
        B = 1e-30;

    C = D/B;

    if(fabs(C)>1e30)
        C=1e30;

    *apg = atan(C);
    if(*apg<0)
        *apg+=pi;
}

void VSAngleSgroup (float vp0, float vs0, float epsilon, float delta, float ap, 
                   float *vsg, float *asg)
/*< VSAngleSgroup : calculate group velocity and angle from phase velocity for
 * given phase angle of qSV-wave in 2D VTI media >*/
{
    float alpha = vs0/vp0;
    float f = 1-alpha*alpha;

    float vsp = VSphase (vp0, vs0, epsilon, delta, ap);
        
    float A = epsilon*sin(2*ap)
                - f*f/4.0/sqrt(4*epsilon*epsilon*pow(sin(ap),4)
                    + f*4.0*sin(ap)*sin(ap)*(2*delta*cos(ap)*cos(ap)-epsilon*cos(2*ap))+f*f)
                * (4*epsilon*sin(2*ap)/f+8*epsilon*epsilon*sin(ap)*sin(ap)*sin(2*ap)/(f*f)
                    - 8*(epsilon-delta)*sin(2*ap)*cos(2*ap)/f);

    float dvda = A*vp0*vp0/vsp/2.0;

    float B, D, C;

    *vsg = vsp*sqrt(1+pow(dvda/vsp,2));

    B=1-tan(ap)*dvda/vsp;
    D=tan(ap)+dvda/vsp;

    if(fabs(B)<1e-30)
        B=1e-30;

    C=D/B;
    if(fabs(C)>1e20)
        C=1e20;

    *asg=atan(C);
    if(*asg<0)
        *asg+=pi;
}

int main(int argc, char* argv[])
{

    int    nx, nz;
    float  dx, dz, time, phaseAngle; 
    float  vp0, vs0, epsilon, delta, theta;
    float  favg;
    int   i,j,d_x,d_z;
    float vpg, apg, vsg, asg, ap, d;
    float *wfp, *wfs, *wf;
    int   it, nt;
    float x, a, dt;
    float amax=0.0;
    int   itmax=0;
    FILE *Fo1, *Fo2, *Fo3;
    int sx;
    int sz; 

    /* setup I/O files */
    Fo1 = fopen("ps.dat","wb"); /* P & SV waves wavefront */
    Fo2 = fopen("pe0.5d0.5.dat","wb"); /* P-wave wavefront */
    Fo3 = fopen("s.dat","wb"); /* SV-wave wavefront */

/*        */ nx=400;
/*        */ nz=400;
/*      m */ dx=5;
/*      m */ dz=5;

/*     vp */ vp0=2000.0;
/*     vs */ vs0=1200.0;

/*epsilon */ epsilon=0.5;
/*  delta */ delta=0.5;
/* degree */ theta=0.0;

/*     Hz */ favg=400.0;

/*        */ nt=601;
/*      s */ dt=0.0005;

/*        */ sx=nx/2;
/*        */ sz=nz/2; 


    time=dt*nt;
    theta *= pi/180.0;
    phaseAngle = 0.05;


    fprintf(stderr,"time= %f \n",time);
        
    for(it=0;it<nt;it++) {

       a = WltRicker(it*dt,favg, 1.0);
	    if(fabs(a)>amax){
	        amax=fabs(a);
	        itmax=it;
	    } 
    }
    time -= itmax*dt;
    fprintf(stderr,"time= %f \n",time);

    wfp = (float*)malloc(nx*nz*sizeof(float));
    wfs = (float*)malloc(nx*nz*sizeof(float));
    wf = (float*)malloc(nx*nz*sizeof(float));

    zero1float(wfp,nz*nx);
    zero1float(wfs,nz*nx);
    zero1float(wf,nz*nx);

    for(i=0;i<nx*nz;i++) {

	        wf[i] = 0.0;
	        wfp[i] = 1;
	        wfs[i] = -1;
    }


    for(i=0;i<360.0/phaseAngle;i++) {

        ap = i*phaseAngle*pi/180.0;

        //vpg = VPgroup (vp0, vs0, epsilon, delta, ap);
        //vsg = VSgroup (vp0, vs0, epsilon, delta, ap);
        //apg = AnglePgroup (vp0, vs0, epsilon, delta, ap);
        //asg = AngleSgroup (vp0, vs0, epsilon, delta, ap);

        VPAnglePgroup (vp0, vs0, epsilon, delta, ap, &vpg, &apg);
        VSAngleSgroup (vp0, vs0, epsilon, delta, ap, &vsg, &asg);

        if(i >= 180.0/phaseAngle) {
            apg += pi;
            asg += pi;
        }

        d=vpg*time;
        d_x=(int)(d*sin(theta+apg)/dx)+sx;
        d_z=-(int)(d*cos(theta+apg)/dz)+sz;
	
        if(d_x>=0 && d_x<nx && d_z>=0 && d_z<nz)
            wfp[d_x*nz+d_z]=-1;

        d=vsg*time;
        d_x=(int)(d*sin(theta+asg)/dx)+sx;
        d_z=-(int)(d*cos(theta+asg)/dz)+sz;
	
        if(d_x>=0 && d_x<nx && d_z>=0 && d_z<nz)
            wfs[d_x*nz+d_z]=1;
    }

    for(i=0;i<nx*nz;i++) {

        wf[i] = wfp[i] + wfs[i];
        wfp[i] = (1-wfp[i])/2;
        wfs[i] = (wfs[i]+1)/2;
        fwrite(&wf[i],sizeof(float), 1, Fo1);
        fwrite(&wfp[i],sizeof(float), 1, Fo2);
        fwrite(&wfs[i],sizeof(float), 1, Fo3);
    }

    free(wfp);
    free(wfs);
    free(wf);
}

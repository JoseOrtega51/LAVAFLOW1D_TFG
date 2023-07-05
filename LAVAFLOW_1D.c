/* Dambreak simulator
   Computational Hydraulics Group
   Codigo base: Pilar García Navarro para Modelos y simulacion de flujos e instalaciones

   José Ortega Moya TFG Física
  =========================================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib/DEFINE.h"
#include "lib/functions.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))		// Macro to find the minimum of two numbers

// File pointers
FILE *dataX_t1;
FILE *dataX_t2;
FILE *time_out;

main()
{

    int i,ic,iter=0,Nout=0;
    double *h,*h0,*A,*A0,*Q,*Q0,*E,*E0,*T,*T0,*r,*r0,*TD,*TD0, *v_rho, *v_rho0;					//  variables
    double *c,*c0,*v,*v0;						// Wave speed
    double *Fr,*Fr0;						// Froude number
    double *DA;						// Variation Conserved variable  A in Dt
    double *DQ;             //Variation Conserved variable  Q in Dt
    double *DE;             //Variation Conserved variable  E in Dt
    double *Sf, *St;             //TErminos fuente
    double *Rh;				// Radio hidraulico

    double lambda1,lambda2,lambda3;				// Eigenvalues
    double e1[3],e2[3],e3[3];					// Eigenvectors
    double alpha1,alpha2,alpha3;		                  	// Wave-strenghts
    double lb1l, lb1r, lb2l, lb2r,lb3l,lb3r, lambda1plus, lambda1minus, lambda2plus, lambda2minus,lambda3plus,lambda3minus;     //Entropy correction

    double t=0.0,dt;
    double dtmin;

    double uPROMEDIO,HPROMEDIO,cPROMEDIO,mPROMEDIO,deltaA,deltaQ,deltaE;    //Variables promediadas Roe

    double L,x;

    double beta=0.;
    double *zb,*delta_zb, *zs;

    double APROMEDIO=0.,TDPROMEDIO=0.,rPROMEDIO=0.;

    double tLastPrint=printTime;        //Control print snapshot

    double tau_b=0.0;

    double A1,A2,B1,B2, n1,n2;

    char name[100];

    int check=0;
    double Ts=0.;


    printf("\n   UPWIND DAMBREAK SIMULATOR");
    printf("\n   Version 1.0");
    printf("\n   --------------------------------------------------");

    i=system ("rm outputFiles/*.out");	// Remove old files
    i=system ("rm OutputTransient/*.txt");	// Remove old files

    dataX_t1=fopen("outputFiles/datosX_t1.out","w");
    dataX_t2=fopen("outputFiles/datosX_t2.out","w");
    time_out=fopen("print_times.txt","w");


    printf("\n>> Simulation setup loaded");
    printf("\n   Simulation time: %.5lf s",simTime);
    printf("\n   CFL: %.2lf",CFL);
    printf("\n   Number of cells: %d",nCells);
    printf("\n   dx: %.2lf m",dx);
    printf("\n   Pipe length: %.2lf m",(nCells-1)*dx);

// Memory allocation
    h = (double *)malloc(nCells*sizeof(double));
    h0 = (double *)malloc(nCells*sizeof(double));
    v = (double *)malloc(nCells*sizeof(double));
    v0 = (double *)malloc(nCells*sizeof(double));
    A = (double *)malloc(nCells*sizeof(double));
    A0 = (double *)malloc(nCells*sizeof(double));
    Q = (double *)malloc(nCells*sizeof(double));
    Q0 = (double *)malloc(nCells*sizeof(double));
    c = (double *)malloc(nCells*sizeof(double));
    c0 = (double *)malloc(nCells*sizeof(double));
    DA = (double *)malloc(nCells*sizeof(double));
    DQ = (double *)malloc(nCells*sizeof(double));
    Fr = (double *)malloc(nCells*sizeof(double));
    Fr0 = (double *)malloc(nCells*sizeof(double));
    Sf = (double *)malloc(nCells*sizeof(double));
    St = (double *)malloc(nCells*sizeof(double));
    Rh = (double *)malloc(nCells*sizeof(double));
    zb = (double *)malloc(nCells*sizeof(double));
    delta_zb = (double *)malloc(nCells*sizeof(double));
    zs = (double *)malloc(nCells*sizeof(double));
    E = (double *)malloc(nCells*sizeof(double));
    E0 = (double *)malloc(nCells*sizeof(double));
    T = (double *)malloc(nCells*sizeof(double));
    T0 = (double *)malloc(nCells*sizeof(double));
    r = (double *)malloc(nCells*sizeof(double));
    r0 = (double *)malloc(nCells*sizeof(double));
    TD = (double *)malloc(nCells*sizeof(double));
    TD0 = (double *)malloc(nCells*sizeof(double));
    DE = (double *)malloc(nCells*sizeof(double));
    v_rho0 = (double *)malloc(nCells*sizeof(double));
    v_rho = (double *)malloc(nCells*sizeof(double));

    printf("\n>> Memory allocation completed");

// Variable initialization
    for(ic=0; ic<nCells; ic++)
    {
        h[ic] = 0.0;
        h0[ic] = 0.0;
        v[ic] = 0.0;
        v0[ic] = 0.0;
        A[ic] = 0.0;
        A0[ic] = 0.0;
        Q[ic] = 0.0;
        Q0[ic] = 0.0;
        c0[ic] = 0.0;
        c[ic] = 0.0;
        Fr[ic] = 0.0;
        Fr0[ic] = 0.0;
        Rh[ic] = 0.0;
        Sf[ic] = 0.0;
        zb[ic] = 0.0;
        delta_zb[ic] = 0.0;
        zs[ic] = 0.0;
        E[ic] = 0.0;
        E0[ic] = 0.0;
        T[ic] = 0.0;
        T0[ic] = 0.0;
        r[ic] = 0.0;
        r0[ic] = 0.0;
        TD[ic] = 0.0;
        TD0[ic] = 0.0;
        DE[ic] = 0.0;
        St[ic] = 0.0;
        v_rho[ic] = 0.0;
        v_rho0[ic] = 0.0;

    }
    printf("\n>> Variable initialization completed");

    FILE *fp = fopen("DATA_FIT_BINGHAM.txt", "r");

    if (fp == NULL)
    {
        printf("Error: could not open file");
        return 1;
    }

    fscanf(fp,"%lf %lf %lf %lf %lf %lf",&A1,&B1,&n1,&A2,&B2,&n2);
    fclose(fp);

// Initial conditions
    L=(nCells-1)*dx;

    for(ic=0; ic<nCells; ic++)
    {

        x=ic*dx;

// CONDICIONES INIICIALES POR TRAMOS.
#if TEST_CASE==1
    if(x<305){
        zs[ic]=30.5;
        Q0[ic]=0.;
        T0[ic]=500.;
    }else{
        zs[ic]=0.;
        Q0[ic]=0.;
        T0[ic]=500.;
    }
    zb[ic]=0.;
#endif // TEST_CASE

#if TEST_CASE==2
        zb[ic]=100.-100./L*x;
        zs[ic]=zb[ic];
        Q0[ic]=0.;
        if(x<100.){
            zs[ic]=zb[ic]+0.;
            //T0[ic]=1500.;
        }
        zs[0]=zb[0]+3.;
        T[0]=1500;
        Q[0]=2.*(1.+K_rho/rho0*(1500.-To));
#endif // TEST_CASE
/**
        if(x<1000.)
        {
            zs[ic]=0.;
            Q0[ic]=0.;
            T0[ic]=323.15;
        }
        else if(x>=1000.&& x<=2000.)
        {
            zs[ic]=0.;
            Q0[ic]=0.0;
            T0[ic]=1000.;
        }else{
            zs[ic]=0.;
            Q0[ic]=0.;
            T0[ic]=1000.;
        }
zs[0]=1.;**/
/**
        if(x>=2.*nCells*dx/5. && x<=3.*nCells*dx/5.)
        {
            zb[ic]=0.;
            //zs[ic]=3.;
        }


/**
        zb[ic]=100.-100./L*x;
        zs[ic]=zb[ic];
        Q0[ic]=0.;
        if(x<500.){
            zs[ic]=zb[ic]+10.;
            T0[ic]=1500.;
        }
/**

        if(x<2000.){
            zs[ic]=zb[ic]+1.;
            T0[ic]=400.;
        }else if(x>=2000.){
            if(!check){
               zs[ic]=zb[ic]+0.1;
            }else{
                zs[ic]=zs[ic-1];
            }
            zb[ic]=0.;
            check=1;
            T0[ic]=200.;
        }
**/

//Inicializacion de las variables
        h0[ic]=zs[ic]-zb[ic];
        if(h0[ic]<TOL)
        {
            v0[ic]=0.0;
            c0[ic]=0.0;			// Wave speed
            Fr0[ic]=0.0;
            TD0[ic]=0.0;
            r0[ic]=0.0;
            E0[ic]=0.0;
            A0[ic]=0.0;
            v_rho0[ic] = 0.0;
        }
        else
        {

            TD0[ic]=K_rho/rho0*(T0[ic]-To);
            r0[ic]=1.+TD0[ic];
            E0[ic]=h0[ic]*TD0[ic];
            c0[ic]=sqrt(0.5*g*h0[ic]*(1+r0[ic]-TD0[ic]));			// Wave speed
            Fr0[ic]=v0[ic]/c0[ic];
            A0[ic]=B*(h0[ic])*r0[ic];
            v0[ic]=Q0[ic]/A0[ic];
            v_rho0[ic] = rho0+K_rho*(T0[ic]-To);

        }


        //Friccion
        if(A0[ic]<TOL)
        {
            Sf[ic]=0.0;
            Rh[ic]=0.0;
            St[ic]=0.0;
        }
        else
        {
            Rh[ic]=h0[ic]*B/(B+2*h0[ic]);
            if(viscous_model==1)
            {
                tau_b=Tb_Newton(v0[ic],h0[ic]);
            }
            else if(viscous_model==2)
            {
                tau_b=Tb_Bing(v0[ic],h0[ic]);
            }
            else if(viscous_model==3)
            {
                tau_b=Tb_Bing_T_mean(v0[ic],h0[ic],T0[ic]);
            }
            else if(viscous_model==4)
            {
                tau_b=Tb_Bing_profile(v0[ic],h0[ic],T0[ic],A1,B1,n1,A2,B2,n2);
            }

            Sf[ic]=n*n*v0[ic]*fabs(v0[ic])/pow(Rh[ic],4./3.);// Manning
            Sf[ic]=Sf[ic]+tau_b/rho0/g/A0[ic]*B; //Modelo viscoso
            /**
            if(h[ic]<TOL){
                St[ic]=0.0;
            }else{St[ic]=Qh(T[ic],T[ic])/rho0/A[ic]/Cp;}
            //St[ic]=Qh(T0[ic])/rho0/A0[ic]/Cp;
            **/

            #if viscous_model==4
                    if(PERFIL==1){
                        Ts=4*T[ic]-3*1500;
                    }else if(PERFIL==2){
                        if(T0[ic]<1499.){
                            Ts=1548*pow((1501-T0[ic]),-0.25);
                        }else Ts=1500.;
                    }
                    #if CONDUCCION_TRAMOS==1
                        if(h[ic]>HMIN){
                            St[ic]=(Qh(Ts,T0[ic])-k_cond*4/h0[ic]*(Ts-T_air))/rho0/A0[ic]/Cp;
                        }else{
                            St[ic]=(Qh(Ts,T0[ic])-k_cond*4/HMIN*(Ts-T_air))/rho0/A0[ic]/Cp;
                        }
                    #else
                        St[ic]=(Qh(Ts,T0[ic]))/rho0/A0[ic]/Cp;
                    #endif // CONDUCCION_TRAMOS
                #else
                    St[ic]=Qh(T0[ic],T0[ic])/rho0/A0[ic]/Cp;
                #endif // viscous_model
        }
    }
    for(ic=0; ic<nCells; ic++)
    {
        T[ic]=TD0[ic]*rho0/K_rho+To;
        fprintf(dataX_t1,"%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n",ic*dx,h0[ic],v0[ic],zb[ic],T[ic]);
    }

    //Delta zb computation
    for(ic=0; ic<nCells-1; ic++)
    {
        delta_zb[ic]=zb[ic+1]-zb[ic];
        //printf("%lf\n",delta_zb[ic]);
    }
    printf("\n>> Initial conditions loaded");

//find dt based on CFL condition
    dt=dx/(fabs(v0[0])+c0[0]);
    for(ic=1; ic<nCells; ic++)
    {

        dtmin=dx/(fabs(v0[ic])+c0[ic]);
        dt=MIN(dt,dtmin);

    }
    dt = dt*CFL;						// Time step computation


// TIME LOOP ///////////////////
    while(t<=simTime)
    {

        for(ic=0; ic<nCells; ic++)
        {
            DA[ic]=0.0;
            DQ[ic]=0.0;
            DE[ic]=0.0;
        }
        for(ic=0; ic<nCells-1; ic++)
        {

//Roe averages

            if (h0[ic]<TOL && h0[ic+1]<TOL)
            {
                uPROMEDIO = 0.0;
                cPROMEDIO = 0.0;
                rPROMEDIO = 0.0;
                mPROMEDIO = 0.0;
                HPROMEDIO = 0.0;
                TDPROMEDIO = 0.0;
            }
            else
            {
                rPROMEDIO = (A0[ic]+A0[ic+1])/(h0[ic]+h0[ic+1])/B;
                /**
                if(A0[ic]<TOL || A0[ic+1]<TOL){
                    TDPROMEDIO=0.0;
                }else{
                    TDPROMEDIO = rPROMEDIO*(TD0[ic]*h0[ic]*sqrt(A0[ic+1])+TD0[ic+1]*h0[ic+1]*sqrt(A0[ic]))/(A0[ic]/B*sqrt(A0[ic+1])+A0[ic+1]/B*sqrt(A0[ic]));
                }**/
                //TDPROMEDIO = rPROMEDIO*(E0[ic]*sqrt(A0[ic+1])+E0[ic+1]*sqrt(A0[ic]))/(A0[ic]/B*sqrt(A0[ic+1])+A0[ic+1]/B*sqrt(A0[ic]));
                TDPROMEDIO = rPROMEDIO*(E0[ic]*sqrt(A0[ic])+E0[ic+1]*sqrt(A0[ic+1]))/(A0[ic+1]/B*sqrt(A0[ic+1])+A0[ic]/B*sqrt(A0[ic]));
                mPROMEDIO = (1.+rPROMEDIO-TDPROMEDIO);
                HPROMEDIO = 0.5*(h0[ic+1]+h0[ic]);
                uPROMEDIO  = (sqrt(A0[ic])*v0[ic]+sqrt(A0[ic+1])*v0[ic+1])/(sqrt(A0[ic])+sqrt(A0[ic+1]));
                cPROMEDIO = sqrt(0.5*g*HPROMEDIO*mPROMEDIO);
            }
               // printf("%lf\n",TDPROMEDIO);
            //printf("%lf\t %lf\n",uPROMEDIO,TDPROMEDIO);
// Eigenvalues
            lb1l=v0[ic]-c0[ic];
            lb1r=v0[ic+1]-c0[ic+1];
            lambda1 = uPROMEDIO-cPROMEDIO;

            lb2l=v0[ic]+c0[ic];
            lb2r=v0[ic+1]+c0[ic+1];
            lambda2 = uPROMEDIO+cPROMEDIO;

            lb3l=v0[ic];
            lb3r=v0[ic+1];
            lambda3 = uPROMEDIO;

// Eigenvectors

            e1[0] =1 ;
            e1[1] = lambda1;

            if(rPROMEDIO<TOL)
            {
                e1[2]=0.0;
                e2[2]=0.0;
            }
            else
            {
                e1[2] = TDPROMEDIO/rPROMEDIO;
                e2[2] = TDPROMEDIO/rPROMEDIO;
            }
            e2[0] = 1;
            e2[1] = lambda2;

            e3[0] = rPROMEDIO;
            e3[1] = rPROMEDIO*uPROMEDIO;
            e3[2] = 1+rPROMEDIO;

// Spatial jump in the conserved variables

            deltaA = A0[ic+1]-A0[ic];
            deltaQ = Q0[ic+1]-Q0[ic];
            deltaE = E0[ic+1]-E0[ic];


// Coefficients
            alpha1 = (deltaA/B*((1+rPROMEDIO)*cPROMEDIO+mPROMEDIO*uPROMEDIO)-deltaQ/B*mPROMEDIO-rPROMEDIO*cPROMEDIO*deltaE)/2./cPROMEDIO/mPROMEDIO;
            alpha2 = (deltaA/B*((1+rPROMEDIO)*cPROMEDIO-mPROMEDIO*uPROMEDIO)+deltaQ/B*mPROMEDIO-rPROMEDIO*cPROMEDIO*deltaE)/2./cPROMEDIO/mPROMEDIO;
            //alpha3 = (deltaE-HPROMEDIO*TDPROMEDIO/rPROMEDIO*deltaA/B)/mPROMEDIO;
            alpha3 = (deltaE-TDPROMEDIO/rPROMEDIO*deltaA/B)/mPROMEDIO;
// Source terms

            beta=g*rPROMEDIO*HPROMEDIO*delta_zb[ic]/2./cPROMEDIO;

// Waves computation

            if (h0[ic]<TOL && h0[ic+1]<TOL)
            {
                alpha1 = 0.0;
                alpha2 = 0.0;
                alpha3 = 0.0;
                beta = 0.0;
            }
            if(beta!=0.)
            {
                //printf("\t beta: %lf \t",beta);
            }


// First wave lambda 1

            if(lb1l*lb1r >= 0.0)  //ondas en la misma direccion
            {

                if(lambda1<=0.)  //hacia la izquierda
                {
                    DA[ic] = DA[ic] + (lambda1*alpha1)*e1[0];
                    DQ[ic] = DQ[ic] + (lambda1*alpha1)*e1[1];
                    DE[ic] = DE[ic] + (lambda1*alpha1)*e1[2];
                    DA[ic]=DA[ic] - beta*e1[0] ;
                    DQ[ic]=DQ[ic] - beta*e1[1]  ;
                    DE[ic]=DE[ic] - beta*e1[2]  ;

                }
                else   //hacia la derecha
                {
                    DA[ic+1] = DA[ic+1] + (lambda1*alpha1)*e1[0];
                    DQ[ic+1] = DQ[ic+1] + (lambda1*alpha1)*e1[1];
                    DE[ic+1] = DE[ic+1] + (lambda1*alpha1)*e1[2];
                    DA[ic+1] = DA[ic+1] - beta*e1[0];
                    DQ[ic+1] = DQ[ic+1] - beta*e1[1];
                    DE[ic+1] = DE[ic+1] - beta*e1[2];
                }

            }
            else
            {
                if(lb1l<0.){
                    lambda1minus=lb1l*(lb1r-lambda1)/(lb1r-lb1l);
                    lambda1plus=lb1r*(lambda1-lb1l)/(lb1r-lb1l);
                }else{
                    lambda1plus=lb1l*(lb1r-lambda1)/(lb1r-lb1l);
                    lambda1minus=lb1r*(lambda1-lb1l)/(lb1r-lb1l);
                }

                DA[ic] = DA[ic] + (lambda1minus*alpha1)*e1[0];
                DQ[ic] = DQ[ic] + (lambda1minus*alpha1)*e1[1];
                DE[ic] = DE[ic] + (lambda1minus*alpha1)*e1[2];
                DA[ic] = DA[ic] - beta*e1[0];
                DQ[ic] = DQ[ic] - beta*e1[1];
                DE[ic] = DE[ic] - beta*e1[2];

                DA[ic+1] = DA[ic+1] + (lambda1plus*alpha1)*e1[0];
                DQ[ic+1] = DQ[ic+1] + (lambda1plus*alpha1)*e1[1];
                DE[ic+1] = DE[ic+1] + (lambda1plus*alpha1)*e1[2];
            }
// Second wave lambda 2

            if(lb2l*lb2r >= 0.0)
            {
                if(lambda2<0.)
                {
                    DA[ic] = DA[ic] + (lambda2*alpha2)*e2[0];
                    DQ[ic] = DQ[ic] + (lambda2*alpha2)*e2[1];
                    DE[ic] = DE[ic] + (lambda2*alpha2)*e2[2];
                    DA[ic] = DA[ic] + beta*e2[0];
                    DQ[ic] = DQ[ic] + beta*e2[1];
                    DE[ic] = DE[ic] + beta*e2[2];

                }
                else if (lambda2>0.)
                {
                    DA[ic+1] = DA[ic+1] + (lambda2*alpha2)*e2[0];
                    DQ[ic+1] = DQ[ic+1] + (lambda2*alpha2)*e2[1];
                    DE[ic+1] = DE[ic+1] + (lambda2*alpha2)*e2[2];
                    DA[ic+1] = DA[ic+1] + beta*e2[0];
                    DQ[ic+1] = DQ[ic+1] + beta*e2[1];
                    DE[ic+1] = DE[ic+1] + beta*e2[2];
                }
            }
            else
            {
                if(lb2l<0.){
                    lambda2minus=lb2l*(lb2r-lambda2)/(lb2r-lb2l);
                    lambda2plus=lb2r*(lambda2-lb2l)/(lb2r-lb2l);
                }else{
                    lambda2plus=lb2l*(lb2r-lambda2)/(lb2r-lb2l);
                    lambda2minus=lb2r*(lambda2-lb2l)/(lb2r-lb2l);
                }
                DA[ic] = DA[ic] + (lambda2minus*alpha2)*e2[0];
                DQ[ic] = DQ[ic] + (lambda2minus*alpha2)*e2[1];
                DE[ic] = DE[ic] + (lambda2minus*alpha2)*e2[2];
                DA[ic] = DA[ic] + beta*e2[0];
                DQ[ic] = DQ[ic] + beta*e2[1];
                DE[ic] = DE[ic] + beta*e2[2];

                DA[ic+1] = DA[ic+1] + (lambda2plus*alpha2)*e2[0];
                DQ[ic+1] = DQ[ic+1] + (lambda2plus*alpha2)*e2[1];
                DE[ic+1] = DE[ic+1] + (lambda2plus*alpha2)*e2[2];
            }

// Third wave lambda 3

            if(lb3l*lb3r >= 0.0)
            {
                if(lambda3<0.)
                {
                    DA[ic] = DA[ic] + (lambda3*alpha3)*e3[0];
                    DQ[ic] = DQ[ic] + (lambda3*alpha3)*e3[1];
                    DE[ic] = DE[ic] + (lambda3*alpha3)*e3[2];

                }
                else if (lambda3>0.)
                {
                    DA[ic+1] = DA[ic+1] + (lambda3*alpha3)*e3[0];
                    DQ[ic+1] = DQ[ic+1] + (lambda3*alpha3)*e3[1];
                    DE[ic+1] = DE[ic+1] + (lambda3*alpha3)*e3[2];
                }
            }
            else
            {
                if(lb3l<0.){
                    lambda3minus=lb3l*(lb3r-lambda3)/(lb3r-lb3l);
                    lambda3plus=lb3r*(lambda3-lb3l)/(lb3r-lb3l);
                }else{
                    lambda3plus=lb3l*(lb3r-lambda3)/(lb3r-lb3l);
                    lambda3minus=lb3r*(lambda3-lb3l)/(lb3r-lb3l);
                }
                DA[ic] = DA[ic] + (lambda3minus*alpha3)*e3[0];
                DQ[ic] = DQ[ic] + (lambda3minus*alpha3)*e3[1];
                DE[ic] = DE[ic] + (lambda3minus*alpha3)*e3[2];

                DA[ic+1] = DA[ic+1] + (lambda3plus*alpha3)*e3[0];
                DQ[ic+1] = DQ[ic+1] + (lambda3plus*alpha3)*e3[1];
                DE[ic+1] = DE[ic+1] + (lambda3plus*alpha3)*e3[2];
            }
        }



// Cell sweep updating
        for(ic=0; ic<nCells; ic++)
        {

// Using waves

            A[ic] = A0[ic] - DA[ic]*dt/dx;


            Q[ic] = Q0[ic] - DQ[ic]*dt/dx-g*A0[ic]*Sf[ic]*dt;


            E[ic] = E0[ic] - DE[ic]*dt/dx;
            //printf("%lf\n",DE[ic]);
            if(Q[ic]*Q0[ic]<0.0)
            {
                //Q[ic] = Q0[ic] - DQ[ic]*dt/dx;
                Q[ic]=0.;
                //printf("\n Q=0 en %lf",ic*dx);
            }
            //if(Q[ic]!=0.0)printf("\n %lf\t%lf\t%d",Q[ic],Sf[ic],ic);
            //ALtura de fluido (Paso intermedio)
            h[ic] = A[ic]/B-E[ic];

            if(h[ic]<TOL)
            {
                v[ic]=0.0;
                Fr[ic]=0.0;
                c[ic]=0.0;
                TD[ic]=0.0;
                r[ic]=0.0;
            }
            else
            {

                v[ic]=Q[ic]/A[ic];      //Velocidad del fluido

                if (fabs(v[ic])<0.00001)
                    v[ic]=0.0;
                if (v[ic]<0.00&&h[ic]<delta_zb[ic]) {
                        Q[ic]=0.0;
                    v[ic]=0.0;
                }


                TD[ic]=E[ic]/h[ic];     //Calculo de la temperatura adimensional
                r[ic]=A[ic]/h[ic]/B;         //Calculo densidad normalizada r
                c[ic]=sqrt(0.5*g*h[ic]*(1+r[ic]-TD[ic]));
                Fr[ic]=v[ic]/c[ic];
            }
            T[ic]=rho0/K_rho*TD[ic]+To; //Temperatura provisional

            T[ic]= T[ic]+St[ic]*dt; //Actualza temperatura con flujo de calor
            //printf("%lf\n",St[ic]*dt);
            if(T[ic]<To){
                T[ic]=To;
            }
            if(h[ic]<TOL){
                v_rho[ic]=0.0;
                h[ic]=0.0;
                v[ic]=0.0;
            }else{
                v_rho[ic]=rho0+K_rho*(T[ic]-To); //Nueva densidad
                h[ic]= rho0*A[ic]/v_rho[ic];    //Nueva altura
                if(h[ic]<TOL){
                    T[ic]=To;
                    h[ic]=0.0;
                    v_rho[ic]=0.0;
                    v[ic]=0.0;
                }
            }


        }

// UPSTREAM BC

        //Q[0]=0.0;
        //A[0] = A0[0] - DA[0]*dt/dx-g*A0[0]*Sf[0]*dt;
        A[0] = A0[0] - DA[0]*dt/dx;
        Q[0]=2.*(1.+K_rho/rho0*(1500.-To));
        if(v[0]<TOL){
            v[0]=1.;
        }
        T[0]=1500.;

        if(A[0]<TOL)
        {
            v[0]=0.0;
            Fr[0]=0.0;
        }
        else
        {
            v[0]=Q[0]/A[0];
            Fr[0]=v[0]/c[0];
        }

// DOWNSTREAM BC

        Q[nCells-1]=0.0;
        A[nCells-1] = A0[nCells-1] - DA[nCells-1]*dt/dx-g*A0[nCells-1]*Sf[nCells-1]*dt;

        if(A[nCells-1]<TOL)
        {
            v[nCells-1]=0.0;
            Fr[nCells-1]=0.0;
        }
        else
        {
            v[nCells-1]=Q[nCells-1]/A[nCells-1];
            Fr[nCells-1]=v[nCells-1]/c[nCells-1];
        }


        iter++;



// Updating time

        tLastPrint=tLastPrint+dt;
        t=t+dt;

        printf("\n%lf  %lf", dt, t);
        #if FIX_H==1
            h[0]=5.;
        #endif // FIX_H

        for(ic=0; ic<nCells; ic++)
        {
            v0[ic]=v[ic];
            h0[ic]=h[ic];


            T0[ic]=T[ic];
            TD0[ic]=K_rho/rho0*(T[ic]-To);
            r0[ic]=1+TD0[ic];

            A0[ic]=h0[ic]*r0[ic];
            Q0[ic]=h0[ic]*r0[ic]*v0[ic];
            E0[ic]=h0[ic]*TD0[ic];
            ;

            v_rho0[ic]=v_rho[ic];
            if(A0[ic]<0.0)
            {
                printf("ERROR: Altura negativa");
            }

            c0[ic]=c[ic];
            Fr0[ic]=Fr[ic];



            if(A0[ic]<TOL)
            {
                Rh[ic]=0.0;
                Sf[ic]=0.0;
                St[ic]=0.0;
            }
            else
            {
                Rh[ic]=h0[ic]*B/(B+2.*h0[ic]);
                if(viscous_model==1)
                {
                    tau_b=Tb_Newton(v0[ic],h0[ic]);
                }
                else if(viscous_model==2)
                {
                    tau_b=Tb_Bing(v0[ic],h0[ic]);
                }
                else if(viscous_model==3)
                {
                    tau_b=Tb_Bing_T_mean(v0[ic],h0[ic],T0[ic]);
                }
                else if(viscous_model==4)
                {
                    tau_b=Tb_Bing_profile(v0[ic],h0[ic],T0[ic],A1,B1,n1,A2,B2,n2);
                }
                if(tau_b>TOL){
                   // printf("%lf \n",tau_b);
                }

                Sf[ic]=n*n*v0[ic]*fabs(v0[ic])/pow(Rh[ic],4./3.); //Manning
                Sf[ic]=Sf[ic]+tau_b/rho0/g/A0[ic]*B; //Modelo visciso
                //if(Sf[ic]!=0.0) printf("\n %lf\t %d",tau_b,ic);
                #if viscous_model==4
                    if(PERFIL==1){
                        Ts=4*T[ic]-3*1500;
                    }else if(PERFIL==2){
                        if(T0[ic]<1499.){
                            Ts=1548*pow((1501-T0[ic]),-0.25);
                        }else Ts=1500.;
                    }
                    #if CONDUCCION_TRAMOS==1
                        if(h[ic]>HMIN){
                            St[ic]=(Qh(Ts,T0[ic])-k_cond*4/h0[ic]*(Ts-T_air))/rho0/A0[ic]/Cp;
                        }else{
                            St[ic]=(Qh(Ts,T0[ic])-k_cond*4/HMIN*(Ts-T_air))/rho0/A0[ic]/Cp;
                        }
                    #else
                        St[ic]=(Qh(Ts,T0[ic]))/rho0/A0[ic]/Cp;
                    #endif // CONDUCCION_TRAMOS
                #else
                    St[ic]=Qh(T0[ic],T0[ic])/rho0/A0[ic]/Cp;
                #endif // viscous_model
                //printf("%lf \n",St[ic]);
            }
        }

        dt=dx/(fabs(v0[0])+c0[0]);
        for(ic=1; ic<nCells; ic++)
        {

            dtmin=dx/(fabs(v0[ic])+c0[ic]);
            dt=MIN(dt,dtmin);

        }
        // Print to file

        if(tLastPrint>=printTime)
        {
            sprintf(name, "OutputTransient/Flow_%02d.txt", Nout);
            printFlow(name,v,h,T,zb,t,time_out);
            tLastPrint=0.;
            Nout+=1;

        }

        dt = dt*CFL;						// Time step computation
//	getchar();
    }
// END TIME LOOP ///////////////

    for(ic=0; ic<nCells; ic++)
    {
        fprintf(dataX_t2,"%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n",ic*dx,A[ic],Q[ic],zb[ic],TD[ic]);
    }

    fclose(dataX_t1);
    fclose(dataX_t2);

// Memory deallocation
    free(A);
    free(A0);
    free(v);
    free(v0);
    free(Q);
    free(Q0);
    free(Fr);
    free(Fr0);
    free(DA);
    free(DQ);
    free(Sf);
    free(Rh);
    free(zb);
    free(delta_zb);
    free(zs);
    free(E);
    free(E0);
    free(T);
    free(T0);
    free(r);
    free(r0);
    free(TD0);
    free(TD);
    free(DE);
    printf("\n>> Memory deallocation completed");

    printf("\n\n>> Simulation completed!");
    printf("\n\n");
}

///////////////////////////////////////////////////////


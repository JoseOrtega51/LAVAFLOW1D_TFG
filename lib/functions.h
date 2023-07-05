#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED



#endif // FUNCTIONS_H_INCLUDED



void printFlow(char nombre[], double v[], double h[],double T[], double zb[],double time, FILE *f_time)
{
    FILE *f;
    f=fopen(nombre, "w");
    if (f==NULL)
    {
        printf("No se pudo abrir el fichero");
    }
    else
    {
        for(int ic=0; ic<nCells; ic++)
        {
            fprintf(f,"%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n",ic*dx,h[ic],v[ic],zb[ic],T[ic]);

        }
        fprintf(f_time,"%.6lf\n",time);
    }
    fclose(f);
}


double cubicSolver(double a, double b,double d)
{
    double D0=b*b;
    double D1=2.*b*b*b+27.*a*a*d;
    double radic=D1*D1-4.*D0*D0*D0;
    double sol=0.;

    //printf("%lf\n",C);
    if(radic>=0.)
    {
        double C=pow((D1-sqrt(radic))/2.,1./3.);
        sol=-1./3./a*(b+C+D0/C);
    }
    else
    {
        double p=-b*b/3./a/a;
        double q=(2.*b*b*b/27./a/a/a+d/a);
        double t=0.;
        t=2*sqrt(-p/3.)*cos(1/3.*acos(3.*q/2./p*sqrt(-3./p)));
        sol=t-b/3./a;
    }

    return sol;
}


double Tb_Bing_T_mean(double u, double h,double T)
{
    double tau_y_T = A_tau+B_tau*exp(C_tau*T);
    if(h<HMIN|| fabs(u)<TOL2)
    {
        //return tau_y_T;
        return 0;
    }
    else
    {   double tau_y_T = A_tau+B_tau*exp(C_tau*T);
        double mu_T = A_mu*exp(B_mu/T);
        double a=2.;
        double b=-3.*(tau_y_T+2.*mu_T*fabs(u)/h);
        double d=tau_y_T*tau_y_T*tau_y_T;
        double Tb=0.;

        if(u>0.)
        {
            Tb=cubicSolver(a,b,d);
        }
        else if(u<0.)
        {
            Tb=-cubicSolver(a,b,d);
        }
        if(isnan(Tb)){
            return u/fabs(u)*1.e10;
        }else return Tb;

    }

}
//double Tb_Bing(double mu, double tau_y, double u, double)
double Tb_Bing(double u, double h)
{
    if(h<TOL || fabs(u)<TOL)
    {
        return 0.;
    }
    else
    {
        double a=2.;
        double b=-3.*(tau_y+2.*mu*fabs(u)/h);
        double d=tau_y*tau_y*tau_y;
        double Tb=0.;

        if(u>0.)
        {
            Tb=cubicSolver(a,b,d);
        }
        else if(u<0.)
        {
            Tb=-cubicSolver(a,b,d);
        }
        return Tb;


        //printf("%lf\n",Tb);
        return Tb;
    }

}


double Tb_Newton(double u,double h)
{
    if(h<TOL)
    {
        return 0;
    }
    else
    {
        return 3*mu*u/h;
    }
}

double Tb_Bing_profile(double u, double h,double T,double A1, double B1, double n1,double A2, double B2, double n2){

    if(h<TOL|| fabs(u)<TOL)
    {
        return 0.;
    }
    else
    {
        double Tb=0.;

        Tb=(fabs(u)/h-(A2*pow(T/Thot,3)+B2*pow(T/Thot,2)+n2))/(A1+B1*pow(T/Thot,n1));
        if(u>0.)
        {
            Tb=fabs(Tb);
        }
        else if(u<0.)
        {
            Tb=-fabs(Tb);
        }
        if(isnan(Tb)){
            return u/fabs(u)*1.e10;
        }else return Tb;


}
}

double Qh(double Ts, double T){
    //Returns heat flux through the surface
    return -eps_rad*sigma_rad*(T*T*T*T-T_air*T_air*T_air*T_air)-hc_conv*(Ts-T_air);
//return -(T-300.); //500.
}

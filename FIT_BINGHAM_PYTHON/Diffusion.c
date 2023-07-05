/* Cálculo de un perfil de temperatura para la lava mediante modelo de difusión de temperatura para TFG de grado en física.
    Autor: José Ortega Moya.
    Adaptado de un código para simulación de flujo de Couette de Pilar García Navarro (CTMF, Universidad de Zaragoza)
   ---------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>

#define N         2000       // Número de celdas
#define simTime   2000000.0     // Tiempo de simulación (s)
#define h         1.       // Altura lava (m)
#define D         h*4.       // Longitud dominio (m)
#define Pe        0.9       // Número de Pecklet [2·dt/(dy*dy)]
#define rho       3000      // Densidad (kg/m3)
#define printerDT 200000.        // Intervalo de escritura a fichero
#define sigma   5.67e-8     // CTE stefan-boltzmann
#define emiss   0.74        // emissivity
#define hc      50.          // CTE conveccion
#define Tair    300.
#define T0      1500.
#define TOL_T   1.
#define Temp_out 4.
FILE *data;
FILE *temp_file;

int main(){

	double dy,dt,t;
	double u[N],y[N], u0[N];
	int i,step,PrintCounter=0;
	char filename[1024];
	double U=1500.; //Temperatura inicial
	double alpha=1./rho/1100;
	double Tmean, T_old;
	int j;

	dy = D/N;
	dt = Pe*dy*dy/(2.0*alpha);

	t=0.0;
	step=0;
    T_old=T0;
// Condiciones iniciales:
	for(i=0;i<N/4.;i++){
		u[i]=U;
		u0[i]=U;
	}
    for(i=round(N/4.);i<N;i++){
		u[i]=Tair;
		u0[i]=Tair;
	}
	for(i=0;i<N;i++){
		y[i]=i*dy;
	}

// Condiciones de contorno:
	u[0]= U;
	u0[0]=U;
	u[N-1]= Tair;
	u0[N-1]=Tair;
    temp_file=fopen("out_T.txt","w");
	data=fopen("outputData.dat","w");
//CAlculo temperatura media
        Tmean=0.0;
        for(j=0;j*dy<h;j++){
            Tmean+=u[j]*dy;
        }
        Tmean=Tmean/h;
        printf("Tmean=%lf \n",Tmean);

    sprintf(filename,"OutputFiles/outputData%d.dat",PrintCounter);	// Write the numbering in the file name with the corresponding extension.
    data=fopen(filename,"w");
    for(i=0;i<N;i++){
        fprintf(data,"%lf %.6e\n",y[i],u[i]);
    }
    fprintf(temp_file,"%lf \n",Tmean);
    T_old=Tmean;
    PrintCounter++;
// Evolución temporal:
	while(t<=simTime){
		for(i=1;i<N-1;i++){
			u[i] =  u0[i] + alpha*dt/dy/dy*(u0[i-1]-2.*u0[i]+u0[i+1]);
		}

        if(u[0]<Tair+TOL_T){
            u[0]=Tair;
        }else{
            u[0]=u[1]-sigma*emiss*(u0[0]*u0[0]*u0[0]*u0[0]-Tair*Tair*Tair*Tair)*dy-hc*(u[0]-Tair)*dy;
        }

        u[N-2]=u[N-3];
		for(i=0;i<N-1;i++){
			u0[i] = u[i];
		}
//Calculo temperatura media
        Tmean=0.0;
        for(j=0;j*dy<h;j++){
            Tmean+=u[j]*dy;
        }
        Tmean=Tmean/h;
// Imprimimos a un archivo para visualizar en GNUplot:
		if( T_old-Tmean>Temp_out ){

			sprintf(filename,"OutputFiles/outputData%d.dat",PrintCounter);	// Write the numbering in the file name with the corresponding extension.
			data=fopen(filename,"w");

			for(i=0;i<N;i++){
				fprintf(data,"%lf %.6e\n",y[i],u[i]);
			}
			fprintf(temp_file,"%lf \n",Tmean);
			T_old=Tmean;
			PrintCounter++;
		}

        printf("Tmean=%lf \n",Tmean);
// Actualizamos el tiempo
		t=t+dt;

	}
	fclose(data);

}

#ifndef DEFINE_H_INCLUDED
#define DEFINE_H_INCLUDED

#endif // DEFINE_H_INCLUDED

#define FLUIDO 1 // 1- LAVA, 2-FLUIDO EJEMPLO

#define g 9.81
#define TOL 0.001
#define TOL2 0.00001
#define HMIN 0.001      //REDUCIR
#define nCells 3000 //600   // Numero de celdas
#define dx 1. //.5       // spaciado
#define CFL 0.9
#define B 1.       // Anchura canal rectangular: DEJAR SIEMPRE EN 1 (NO SE HA PROBADO CON OTRO VALOR)
#define Q0L 0.0    // caudal inicial izquierda presa
#define Q0R 0.0    // caudal inicial derecha presa
#define h0L 5.0    // calado inicial izquierda presa
#define h0R 1.0    // calado inicial derecha presa
#define simTime 3.6e3//.5//1000.
#define printTime 50// 10.
//#define rho 1835. //Densidad constante (kg/m^3)
#define sigma_rad 6.67e-8

//CONFIGURACION
//TRANSFERENCIA DE CALOR
#define TRANSF_CALOR 1 //1-Si, 2-No
#define CONDUCCION_TRAMOS 1

//CASO TEST PREDEFINIDO
#define TEST_CASE 2 //1-BM1, 2- RAMPA

//MODELO DE FRICCION
#define viscous_model 4 // 1-Newton 2-Bingham 3-Bing T promedio 4-Bing T perfil
#define PERFIL 1 //1- lineal a tramos, 2-calculado por conducción
#define n 0.00  //Coeficiente de manning

//ALTURA FIJADA A LA IZQUIRDA?
#define FIX_H 0

//Propiedades Bingham


#if FLUIDO==1
    #define mu 100. //Viscosidad constante (Pa·s)
    #define tau_y 1500. //Esfuerzo umbral (Pa)
#elif FLUIDO==2
    #define mu 1. //Viscosidad constante (Pa·s)
    #define tau_y 50. //Esfuerzo umbral (Pa)
#elif FLUIDO==3
    #define mu 1000. //Viscosidad constante (Pa·s)
    #define tau_y 1000. //Esfuerzo umbral (Pa)
#endif // FLUIDO


//Dependedncia de la densidad con respecto a la temperatura

#if FLUIDO==1
    #define K_rho -0.1
    #define rho0 3000.
    #define To 300.
#elif FLUIDO==2
    #define K_rho -0.66
    #define rho0 860.84
    #define To 293.15
#elif FLUIDO==3
    #define K_rho -0.1
    #define rho0 1835
    #define To 300
#endif // FLUIDO

//DEPENDENCIA BINGHAM CON LA TEMPERTAURA

#if FLUIDO==1
//LAVA:
    #define A_tau 0.        //Pa
    #define B_tau 5.6e6   //Pa
    #define C_tau -0.0058   //K^-1
    #define A_mu 1.77    //Pa s
    #define B_mu 9500.    //K

    #define Thot 1500.

//Specific heat
    #define Cp 1050.

//HEAT TRANSFERENCE
    #define T_air 300.
    #define eps_rad 0.74
    #define sigma_rad 6.67e-8
    #define hc_conv 50
    #define k_cond 1.
#elif FLUIDO==2
//FLUIDO EJEMPLO
    #define A_tau 0.        //Pa
    #define B_tau 1000.    //Pa
    #define C_tau -0.02   //K^-1
    #define A_mu 0.0067379469990      //Pa s
    #define B_mu 2000.   //K

    #define Thot 323.15

    //Specific heat
    #define Cp 1900.

    //HEAT TRANSFERENCE
    #define T_air 393.15
    #define eps_rad 0.
    #define sigma_rad 6.67e-8
    #define hc_conv 0.

#elif FLUIDO==3

     #define A_tau 0.        //Pa
    #define B_tau  30431.1    //Pa
    #define C_tau -0.00602   //K^-1
    #define A_mu 3.16    //Pa s
    #define B_mu 1726.9   //K

    #define Thot 500.
    #define Cp 1.
     #define T_air 450.
     #define eps_rad 0.74
     #define hc_conv 100.
#endif // FLUIDO

#if TRANSF_CALOR==0
    #define eps_rad 0.
    #define hc_conv 0.
#endif // TRANSF_CALOR



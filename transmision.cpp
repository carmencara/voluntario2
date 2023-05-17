//---------------------------------------------------------------------
//                                                                    |
//          ALGORITMO PARA RESOLVER LA ECUACIÓN DE SCHRÖDINGER        |
//                   PROBLEMA VOLUNTARIO DEL TEMA 3                   |
//                              Objetivos:                            |
//                                                                    |
//   1. Estimar el coeficiente de transmisión K simulando el sistema  |
//      1E3 veces.                                                    |
//   2. Comparar el valor de K y el de p_D(n_D) en t=n_D. Explicar    |
//      la relación observada con argumentos físicos y por qué el     |
//      método de maximización tiene sentido.                         |
//   3. Estimar el error de K.                                        |
//   4. Estudiar la dependencia de K en N (N=500,1000,2000)           |
//   5. Estudiar la dependencia de K en V(x) (λ=0.1,0.3,0.5,1,5,10)   |
//   6. Comparar con resultados teóricos.                             |
//   7. Calcular los valores esperados de observables (posición,      |
//      momento, energía cinética, energía total) con errores.        |
//                                                                    |
//---------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <fstream>     // Para trabajar con ficheros
#include <complex>     // Para trabajar con números complejos
#include "gsl_rng.h"   // Librería para generación de números aleatorios
#include <sys/time.h>

#define N 1000         // Tamaño del retículo espacial
#define PI 3.141593

using namespace std;

// FUNCIONES QUE SE VAN A UTILIZAR
long int SemillaTiempo();
double Random_double(gsl_rng *tau);
complex<double> CalculaPhi(int j, double k0_tilde);
double CalculaNorma2(complex<double> z);

/*---------------------------------------------------------------------
|                           FUNCIÓN PRINCIPAL                         |
---------------------------------------------------------------------*/
int main()
{
    // ------------------ DECLARACIÓN DE VARIABLES --------------------

    int ciclos;                 // Número de oscilaciones completas de la función de onda
    int iteraciones;            // Número de iteraciones que se realizan del algoritmo
    int j,n,k;                  // Contadores
    int mT;                     // Núm. de veces que se detecta a la partícula a la dcha. de la barrera
    int nD;                     // Tiempo para que aparezca un máximo a la dcha. de la barrera
    int t;                      // Tiempo (discretizado)
    int detectar;               // Ejecuta el bucle cuando no se detecta la partícula (detectar=0),
                                // cuando se detecta (detectar=1) el bucle para

    double lambda;              // Cte de proporcionalidad para la energía del fotón incidente
    double s_tilde;             // s/h^2=1/4k0_tilde^2
    double k0_tilde;            // k0*h, con k0*N*h=2*PI*ciclos
    double V[N+1];              // Potencial, es lambda*k0_tilde^2 si j\in[2N/5,3N/5]
    double Aplus, Aminus;       // Para calcular alpha, Aplus=1, Aminus=1
    double norma_phi;           // Norma de la función de onda
    double norma_txt;
    double densidad;            // Densidad de probabilidad: |phi|^2
    double PR;                  // Probabilidad de detectar la partícula a la derecha (right)
    double PL;                  // Probabilidad de detectar la partícula a la izquierda (left)
    double p;                   // Número aleatorio entre 0 y 1
    double K;                   // Coeficiente de transmisión: K=mT/iteraciones
    double K_teorico;           // Coeficiente de transmisión teórico
    double PR_media;            // Probabilidad promedio de hallar la partícula a la derecha
    double varianza;            // Varianza de K
    double error;               // Error estadístico de K = sqrt(varianza/iteraciones)
    double aux;                 // Variable auxiliar
    double h;                   // Paso entre puntos = 2*k0_tilde (para integración numérica)

    complex<double> phi[N+1];   // Función de onda
    complex<double> dphi[N+1];  // Primera derivada de la función de onda
    complex<double> d2phi[N+1]; // Segunda derivada de la función de onda
    complex<double> xi[N+1];    // phi[j,n+1]=xi[j,n]-phi[j,n]
    complex<double> alpha[N];   // xi[j+1]=alpha[j]*xi[j]+beta[j]
    complex<double> beta[N];    // xi[j+1]=alpha[j]*xi[j]+beta[j]
    complex<double> gamma[N];   // gamma[j]=1/(A0[j]+Aplus*alpha[j])
    complex<double> b[N];       // b[j]=2i*phi[j]/s_tilde
    complex<double> A0[N];      // Para calcular alpha, A0[j]=-2+2i/s_tilde-V[j]
    complex<double> expval_x;   // Valor esperado de la posición
    complex<double> error_x;    // Error del valor esperado de la posición
    complex<double> expval_p;   // Valor esperado del momento
    complex<double> error_p;    // Error del valor esperado del momento
    complex<double> expval_T;   // Valor esperado de la energía cinética
    complex<double> error_T;    // Error del valor esperado de la energía cinética
    complex<double> expval_E;   // Valor esperado de la enería total
    complex<double> error_E;    // Error del valor esperado de la energía total

    ofstream fich_potencial;    // Fichero para guardar la forma del potencial
    ofstream fich_K;            // Fichero para guardar todo lo relativo al coef. de transmisión
    ofstream fich_posicion;     // Fichero para guardar el valor esperado de la posición
    ofstream fich_momento;      // Fichero para guardar el valor esperado del momento
    ofstream fich_cinetica;     // Fichero para guardar el valor esperado de la energía cinética
    ofstream fich_energia;      // Fichero para guardar el valor esperado de la energía total


    // ------------------ NÚMEROS ALEATORIOS CON GSL ------------------
    gsl_rng *tau;
    int semilla;
    semilla = SemillaTiempo();
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);


    // ------------------------ INICIALIZACIÓN ------------------------

    fich_potencial.open("potencial.txt");
    fich_K.open("coeficiente_transmision.txt");
    fich_posicion.open("expval_posicion.txt");
    fich_momento.open("expval_momento.txt");
    fich_cinetica.open("expval_cinetica.txt");
    fich_energia.open("expval_energia.txt");

    lambda = 0.9;
    ciclos = 50;                // Restringido a 1,...,N/4
    k0_tilde = 2*PI*ciclos/N;
    s_tilde = 1/(4*k0_tilde*k0_tilde);
    nD = 500;
    h = 2*k0_tilde;

    // Inicialización del potencial
    for(j=0; j<N-1; j++)
    {
        if (j>(2*N/5)&&j<(3*N/5)) V[j] = lambda*k0_tilde*k0_tilde;
        else V[j] = 0.0;
        fich_potencial << j << " " << V[j] << endl;
    }

    // ----------------------- CÁLCULO DE ALPHA -----------------------

    // Inicializar alpha
    alpha[0] = (0.0, 0.0);
    alpha[N] = (0.0, 0.0);

    // Calcular las A
    Aplus = 1.0;
    Aminus = 1.0;
    for(j=1; j<N; j++) A0[j] = complex<double>(-2.0-V[j], 2.0/s_tilde);

    // Calcular alpha y gamma empezando en N-1
    gamma[N-1] = 1.0/A0[N-1];
    for(j=N-2; j>=0; j--)
    {
        alpha[j] = -Aminus*gamma[j+1];
        gamma[j] = 1.0/(A0[j]+alpha[j]);
    }

    // -------------------------- ALGORITMO ---------------------------

    beta[N] = 0.0;
    iteraciones = 1; // (1000 para los apartados 1-6, 1 para el apartado 7)
    mT = 0.0; // Inicializo el número de veces que se detecta a la partícula a la derecha
    PR_media = 0.0; // Inicializo la probabilidad media

    for(n=1;n<=iteraciones;n++)
    {    
        // -- CÁLCULO DE LA FUNCIÓN DE ONDA INICIAL EN CADA ITERACIÓN --

        phi[0] = (0.0, 0.0);
        phi[N] = (0.0, 0.0);
        norma_phi = 0.0;
        for(j=1; j<N; j++)
        {
            phi[j] = CalculaPhi(j, k0_tilde);
            norma_phi = norma_phi + CalculaNorma2(phi[j]);
        }
        norma_phi = sqrt(norma_phi);
        for(j=1; j<N; j++) phi[j]=phi[j]/norma_phi; // Normalizar

        // ----------------- BÚSQUEDA DE LA PARTÍCULA ------------------

        t = 1; // Inicializo el tiempo
        detectar = 0; // Inicialmente la partícula no se detecta (está a la izquierda)

        while(detectar==0)
        {
            // ACTUALIZO EL TIEMPO
            t++;

            // CÁLCULO DE B: b[j]=4i*phi[j]/s_tilde
            for(j=1;j<N;j++) b[j] = complex<double>(0,4)*phi[j]/s_tilde;

            // CÁLCULO DE BETA: beta[j-1]=gamma[j]*(b[j]-Aplus*beta[j])
            for(j=N-2;j>=0;j--) beta[j] = gamma[j+1]*(b[j+1]-beta[j+1]);

            // CÁLCULO DE XI: xi[j+1]=alpha[j]*xi[j]+beta[j]
            xi[0] = 0;
            xi[N] = 0;
            for(j=1;j<N;j++) xi[j] = alpha[j-1]*xi[j-1]+beta[j-1];

            // CÁLCULO DE PHI: phi[j,n+1]=xi[j,n]-phi[j,n]
            norma_phi = 0;
            for(k=1;k<N;k++) phi[k] = xi[k]-phi[k];

            // BUSCO EL MÁXIMO LOCAL DE PD
            if((t%nD)==0)
            {
                // Probabilidad de detectar la partícula a la derecha
                PR = 0.0;
                for(j=0.8*N;j<N;j++) PR = PR + CalculaNorma2(phi[j]);
                PR_media = PR_media + PR;

                // Genero un número real aleatorio entre 0 y 1 para medir a la derecha
                p = Random_double(tau);

                // Partícula encontrada
                if(p<PR)
                {
                    mT++;         // Número de veces que se encuentra
                    cout << mT << " ";
                    detectar = 1; // Termino el bucle
                }
                // Partícula no encontrada
                else
                {
                    // Reinicio phi y su norma
                    for(j=0.8*N;j<N;j++) phi[j]=0;
                    norma_txt=0;
                    for(j=1;j<N;j++) norma_txt = norma_txt + CalculaNorma2(phi[j]);
                    norma_txt=sqrt(norma_txt);
                    for(j=1;j<N;j++) phi[j]=phi[j]/norma_txt;

                    // Probabilidad de detectar la partícula a la izquierda
                    PL = 0.0;
                    for(j=1;j<=0.2*N;j++) PL = PL + CalculaNorma2(phi[j]);
                    
                    // Genero número real aleatorio entre 0 y 1 para medir a la izquierda
                    p = Random_double(tau);

                    // Partícula encontrada
                    if(p<PL)
                    {
                        detectar = 1; // Termino el bucle
                    }
                    // Partícula no encontrada
                    else
                    {
                        for(j=1;j<=N/5;j++) phi[j]=0;
                        norma_txt=0;
                        for(j=1;j<N;j++) norma_txt=norma_txt + CalculaNorma2(phi[j]);
                        norma_txt=sqrt(norma_txt);
                        for(j=1;j<N;j++) phi[j]=phi[j]/norma_txt;
                    }
                }
            }

            // --------------------- VALORES ESPERADOS ---------------------
            aux = exp(-8.0*(4*(N-1)-N)*(4*(N-1)-N)/(1.0*N*N));
            
            // 1. Calcular la primera y segunda derivada de phi
            dphi[0] = phi[1]/aux;
            dphi[N] = (phi[N]-phi[N-1])/h;
            d2phi[0] = 0;
            d2phi[N] = 0;
            for (j=1;j<N;j++)
            {
                dphi[j] = (phi[j+1]-phi[j-1])/(2*h);
                d2phi[j] = (phi[j+1]-2.0*phi[j]+phi[j-1])/(h*h);
            }

            // 2. Inicializar los valores esperados y sus errores
            expval_x = 0;
            expval_p = 0;
            expval_T = 0;
            expval_E = 0;
            error_x = 0;
            error_p = 0;
            error_T = 0;
            error_E = 0;

            // 3. Calcular valores esperados usando integración numérica y errores
            for(j=0;j<=N;j++)
            {
                expval_x = expval_x + j*h*(conj(phi[j])*phi[j]);
                expval_p = expval_p + (conj(phi[j])*dphi[j]);
                expval_T = expval_T + (conj(phi[j])*d2phi[j]);
                expval_E = expval_E + (-conj(phi[j])*d2phi[j])+complex<double>(V[j],0)*(conj(phi[j])*phi[j])/(h*h);

                error_x = error_x + (j*j*h*h*(conj(phi[j])*phi[j]));
                error_p = error_p + (conj(phi[j])*dphi[j])*(conj(phi[j])*dphi[j]);
                error_T = error_T + (conj(phi[j])*d2phi[j])*(conj(phi[j])*d2phi[j]);
                error_E = error_E + (-(conj(phi[j])*d2phi[j])+complex<double>(V[j],0)*(conj(phi[j])*phi[j])/(h*h))*(-(conj(phi[j])*d2phi[j])+complex<double>(V[j],0)*(conj(phi[j])*phi[j])/(h*h));
            }
            error_x = sqrt(error_x - expval_x*expval_x);
            error_p = sqrt(error_p - expval_p*expval_p);
            error_T = sqrt(error_T - expval_T*expval_T);
            error_E = sqrt(error_E - expval_E*expval_E);

            // 4. Guardar los valores esperados en ficheros en función del tiempo
            fich_posicion << t << " " << real(expval_x) << " " << 0 << " " << real(error_x) << endl;
            fich_momento << t << " " << imag(expval_p) << " " << 0 << " " << real(error_p) << endl;
            fich_cinetica << t << " " << sqrt(CalculaNorma2(expval_T)) << " " << 0 << " " << sqrt(CalculaNorma2((error_T))) << endl;
            fich_energia << t << " "<<sqrt(CalculaNorma2(expval_E)) << " " << 0 << " " << sqrt(CalculaNorma2(error_E)) << endl;
        }
    }

    // -------------- COEFICIENTE DE TRANSMISIÓN TEÓRICO --------------
    if(lambda<1) K_teorico = 4*(1-lambda)/(4*(1-lambda)+pow(lambda*sin(2*PI*ciclos*sqrt(1-lambda)/5),2));
    else K_teorico = 4*(lambda-1)/(4*(lambda-1)+pow(lambda*sinh(2*PI*ciclos*sqrt(lambda-1)/5),2));

    // ------------------ COEFICIENTE DE TRANSMISIÓN ------------------ 
    // Coeficiente de transmisión   
    K = 1.0*mT/iteraciones;
    fich_K << "Valor de K:" << "\t"<< K << endl;
    // Coeficiente de transmisión teórico
    fich_K << "K teórico:" << "\t" << K_teorico << endl;
    fich_K << "Error relativo" << "\t" << 100*(K-K_teorico)/K_teorico << " %" << endl;
    // Probabilidad de detectar a la derecha media (para comparar con K)
    PR_media = PR_media/iteraciones;
    fich_K << "P_D promedio:" << "\t" << PR_media <<endl;
    // Error estadístico del coeficiente de transmisión
    varianza = K*(1-K); // Varianza en una distribución binomial
    error = 2.0*sqrt(varianza/iteraciones);
    fich_K << "Error:" << "\t" << error <<endl;

    // Cerrar los ficheros y finalizar
    fich_potencial.close();
    fich_K.close();
    fich_posicion.close();
    fich_momento.close();
    fich_cinetica.close();
    fich_energia.close();
    return 0;
}


/*---------------------------------------------------------------------
|                            OTRAS FUNCIONES                          |
---------------------------------------------------------------------*/

// Función CalculaPhi: devuelve el valor de la componente de Phi
complex<double> CalculaPhi(int j, double k0_tilde)
{
    double aux;
    complex<double> componente_phi;

    aux = exp(-8.0*(4*j-N)*(4*j-N)/(1.0*N*N));
    componente_phi = complex<double>(cos(k0_tilde*j),sin(k0_tilde*j))*aux;    

    return componente_phi;
}

// Función CalculaNorma2: devuelve el valor de la norma al cuadrado de un
// número complejo double
double CalculaNorma2(complex<double> z)
{
    double norma;
    norma = real(z)*real(z) + imag(z)*imag(z);

    return norma;
}


/*---------------------------------------------------------------------
|                  FUNCIONES PARA NÚMEROS ALEATORIOS                  |
---------------------------------------------------------------------*/

// Función que crea una semilla para generar números aleatorios, basada en el tiempo actual
long int SemillaTiempo()
{
    struct timeval tv;
    gettimeofday(&tv,0);
    return (tv.tv_sec + tv.tv_usec);
}

// Función que genera un número REAL ALEATORIO en un rango [0,1]
double Random_double(gsl_rng *tau)
{
    // Generar un número aleatorio entre 0 y 1
    double numero_aleatorio = gsl_rng_uniform(tau);

    return numero_aleatorio;

}
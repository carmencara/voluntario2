//---------------------------------------------------------------------
//                                                                    |
//          ALGORITMO PARA RESOLVER LA ECUACIÓN DE SCHRÖDINGER        |
//                   PROBLEMA VOLUNTARIO DEL TEMA 3                   |
//                              Objetivos:                            |
//                                                                    |
//   8. Estudiar la dependencia de K con el número n de barreras y    |
//      cuantificar este comportamiento haciendo un ajuste con una    |
//      función de prueba.                                            |
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
complex<double> CalculaPhi(int j, double k0_tilde, double Q);
double CalculaNorma2(complex<double> z);

/*---------------------------------------------------------------------
|                           FUNCIÓN PRINCIPAL                         |
---------------------------------------------------------------------*/
int main()
{
    // ------------------ DECLARACIÓN DE VARIABLES --------------------

    int ciclos;                 // Número de oscilaciones completas de la función de onda
    int Q;                      // Ancho del pozo (ahora no es N)
    int barreras;               // Número de barreras
    int iteraciones;            // Número de iteraciones que se realizan del algoritmo
    int i,j,n,k;                // Contadores
    int mT;                     // Núm. de veces que se detecta a la partícula a la dcha. de la barrera
    int nD;                     // Tiempo para que aparezca un máximo a la dcha. de la barrera
    int t;                      // Tiempo (discretizado)
    int detectar;               // Ejecuta el bucle cuando no se detecta la partícula (detectar=0),
                                // cuando se detecta (detectar=1) el bucle para

    // Defino el número de barreras y el ancho del pozo
    barreras = 10;
    Q = N + 2*N*(barreras-1)/5;

    double lambda;              // Cte de proporcionalidad para la energía del fotón incidente
    double s_tilde;             // s/h^2=1/4k0_tilde^2
    double k0_tilde;            // k0*h, con k0*N*h=2*PI*ciclos
    double V[Q+1];              // Potencial, es lambda*k0_tilde^2 si j\in[2N/5,3N/5]
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

    complex<double> phi[Q+1];   // Función de onda
    complex<double> xi[Q+1];    // phi[j,n+1]=xi[j,n]-phi[j,n]
    complex<double> alpha[Q];   // xi[j+1]=alpha[j]*xi[j]+beta[j]
    complex<double> beta[Q];    // xi[j+1]=alpha[j]*xi[j]+beta[j]
    complex<double> gamma[Q];   // gamma[j]=1/(A0[j]+Aplus*alpha[j])
    complex<double> b[Q];       // b[j]=2i*phi[j]/s_tilde
    complex<double> A0[Q];      // Para calcular alpha, A0[j]=-2+2i/s_tilde-V[j]

    ofstream fich_potencial;    // Fichero para guardar la forma del potencial
    ofstream fich_K;            // Fichero para guardar todo lo relativo al coef. de transmisión


    // ------------------ NÚMEROS ALEATORIOS CON GSL ------------------
    gsl_rng *tau;
    int semilla;
    semilla = SemillaTiempo();
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);


    // ------------------------ INICIALIZACIÓN ------------------------

    fich_potencial.open("b_potencial.txt");
    fich_K.open("b_coeficiente_transmision.txt");

    lambda = 0.5;
    ciclos = Q/4;                // Restringido a 1,...,N/4
    k0_tilde = 2*PI*ciclos/Q;
    s_tilde = 1/(4*k0_tilde*k0_tilde);
    nD = round(2000*(barreras/2.0));
    h = 2*k0_tilde;

    // Inicialización del potencial
    // Zona I) Izquierda de las barreras
    for(j=0;j<=N/5;j++)
    {
        V[j] = 0.0; 
        fich_potencial << j << " " << V[j] << endl; 
    }
    // Zona II) Barreras
    for(j=0;j<2*barreras+1;j++)
    { 
        for(k=0;k<=N/5;k++) // La zona de las barreras sigue midiendo N/5
        {
            if(j%2!=0) // Los tramos impares son las zonas altas de las barreras
            {
                i = k + (j+1)*N*0.2;
                V[i] = lambda*k0_tilde*k0_tilde;
                fich_potencial << i << " " << V[i] << endl;
            }
            else // Los tramos pares son las zonas bajas de las barreras
            {
                i = k + (j+1)*N*0.2;
                V[i] = 0.0; 
                fich_potencial << i << " " << V[i] << endl;   
            }
        }
    }
    // Zona III) Derecha de las barreras
    for(j=(2*barreras+1)*N/5;j<=Q;j++)
    {
        V[j] = 0.0; 
        fich_potencial << j << " " << V[j] << endl; 
    }

    // ----------------------- CÁLCULO DE ALPHA -----------------------

    // Inicializar alpha
    alpha[0] = (0.0, 0.0);
    alpha[Q] = (0.0, 0.0);

    // Calcular las A
    Aplus = 1.0;
    Aminus = 1.0;
    for(j=1; j<Q; j++) A0[j] = complex<double>(-2.0-V[j], 2.0/s_tilde);

    // Calcular alpha y gamma empezando en N-1
    gamma[Q-1] = 1.0/A0[Q-1];
    for(j=Q-2; j>=0; j--)
    {
        alpha[j] = -Aminus*gamma[j+1];
        gamma[j] = 1.0/(A0[j]+alpha[j]);
    }

    // -------------------------- ALGORITMO ---------------------------

    beta[Q] = 0.0;
    iteraciones = 1000;
    mT = 0.0; // Inicializo el número de veces que se detecta a la partícula a la derecha
    PR_media = 0.0; // Inicializo la probabilidad media

    for(n=1;n<=iteraciones;n++)
    {    
        // -- CÁLCULO DE LA FUNCIÓN DE ONDA INICIAL EN CADA ITERACIÓN --

        phi[0] = (0.0, 0.0);
        phi[Q] = (0.0, 0.0);
        norma_phi = 0.0;
        for(j=1; j<Q; j++)
        {
            phi[j] = CalculaPhi(j, k0_tilde, Q);
            norma_phi = norma_phi + CalculaNorma2(phi[j]);
        }
        norma_phi = sqrt(norma_phi);
        for(j=1; j<Q; j++) phi[j]=phi[j]/norma_phi; // Normalizar

        // ----------------- BÚSQUEDA DE LA PARTÍCULA ------------------

        t = 1; // Inicializo el tiempo
        detectar = 0; // Inicialmente la partícula no se detecta (está a la izquierda)

        while(detectar==0)
        {
            // ACTUALIZO EL TIEMPO
            t++;

            // CÁLCULO DE B: b[j]=4i*phi[j]/s_tilde
            for(j=1;j<Q;j++) b[j] = complex<double>(0,4)*phi[j]/s_tilde;

            // CÁLCULO DE BETA: beta[j-1]=gamma[j]*(b[j]-Aplus*beta[j])
            for(j=Q-2;j>=0;j--) beta[j] = gamma[j+1]*(b[j+1]-beta[j+1]);

            // CÁLCULO DE XI: xi[j+1]=alpha[j]*xi[j]+beta[j]
            xi[0] = 0;
            xi[Q] = 0;
            for(j=1;j<Q;j++) xi[j] = alpha[j-1]*xi[j-1]+beta[j-1];

            // CÁLCULO DE PHI: phi[j,n+1]=xi[j,n]-phi[j,n]
            //norma_phi = 0;
            for(k=1;k<Q;k++) phi[k] = xi[k]-phi[k];

            // BUSCO EL MÁXIMO LOCAL DE PD
            if((t%nD)==0)
            {
                // Probabilidad de detectar la partícula a la derecha
                PR = 0.0;
                for(j=(N/5)*2*(barreras+1);j<Q;j++) PR = PR + CalculaNorma2(phi[j]);

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
                    for(j=2*(barreras+1)*N/5;j<Q;j++) phi[j]=0;
                    norma_txt=0;
                    for(j=1;j<Q;j++) norma_txt = norma_txt + CalculaNorma2(phi[j]);
                    norma_txt=sqrt(norma_txt);
                    for(j=1;j<Q;j++) phi[j]=phi[j]/norma_txt;

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
                        for(j=1;j<Q;j++) norma_txt=norma_txt + CalculaNorma2(phi[j]);
                        norma_txt=sqrt(norma_txt);
                        for(j=1;j<Q;j++) phi[j]=phi[j]/norma_txt;
                    }
                }
            }
            if(t>4*nD) detectar=1;
        }
    }

    // ------------------ COEFICIENTE DE TRANSMISIÓN ------------------ 
    // Coeficiente de transmisión   
    K = 1.0*mT/iteraciones;
    fich_K << "Valor de K:" << "\t"<< K << endl;
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
    return 0;
}


/*---------------------------------------------------------------------
|                            OTRAS FUNCIONES                          |
---------------------------------------------------------------------*/

// Función CalculaPhi: devuelve el valor de la componente de Phi
complex<double> CalculaPhi(int j, double k0_tilde, double Q)
{
    double aux;
    complex<double> componente_phi;

    aux = exp(-8.0*(4*j-Q)*(4*j-Q)/(1.0*Q*Q));
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
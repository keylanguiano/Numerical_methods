//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float* bairstow (int n, float *a, float *b, float *db, float e_a, int ni, float r, float s);

float* bairstow (int n, float *a, float *b, float *db, float e_a, int ni, float r, float s)
{
    int i, j, k, l, m, o;

    float R[2], vr[ni], vs[ni], det, dr, ds, raiz[n], ea, ear, eas, *aux, *xr;

    xr = (float*) calloc(n * 2, sizeof(float));

    // Continuidad del programa de acuerdo al numero de raices necesarias
    for (i = 0; i < (n * 2) ; i+=4)
    {
        // Asignacion de valores iniciales para r y s
        vr[0] = r;
        vs[0] = s;

        // Presicion de la raiz de acuerdo al numero maximo de iteraciones
        for(j = 0, ea = 100 ; j < ni; j++)
        {
            b[n] = a[n];
            b[n-1] = a[n-1] + b[n] * vr[j];

            for(k = (n - 2); k > -1; k--)
            {
                b[k] = a[k] + b[k+1] * vr[j] + b[k+2] * vs[j];
            }

            // Termino lineal
            R[1] = b[1];

            // Termino independiente
            R[0] = b[0] - b[1] * vr[j];
        
            db[n] = 0;
            db[n - 1] = b[n];

            for(l = (n - 2); l > -1; l--)
            {
                db[l] = db[l+1] * vr[j] + db[l+2] * vs[j] + b[l+1];
            }

            det = 1 / (db[2] * db[0] - db[1] * db[1]);
            
            dr = det * (db[1] * b[1] - db[2] * b[0]);
            ds = det * (db[1] * b[0] - db[0] * b[1]);
            
            // Valores nuevos para r y s
            vr[j + 1] = vr[j] + dr;
            vs[j + 1] = vs[j] + ds;

            // Calculo de los errorres de aproximacion por r y s
            ear = fabs (dr / vr[j + 1]);
            eas = fabs (ds / vs[j + 1]);

            // Seleccion del error de aproximacion para la iteracion
            if (eas > ea)
            {
                ea = eas;
            }else
            {
                ea = ear;
            }

            printf("\n%d. R[1] = %f\tR[0] = %f (%f)", j+1, R[1], R[0], ea);

            // Comprobacion de la continuidad del programa por error de aproximacion
            if(ea < e_a)
            {
                break;
            }
        }

        printf("\n");

        // Asignacion de las raices segun sus valores reales e imaginarios en el vector segun la paridad
        raiz[i] = vr[j] * vr [j]+ 4 * vs[j];
        
        xr[i] = (vr[j] + (raiz[i] < 0 ? 0 : sqrt(raiz[i]))) / 2;
        xr[i+1] = raiz[i] < 0? sqrt(-raiz[i]) / 2 : 0;

        xr[i+2] = (vr[j] - (raiz[i] < 0 ? 0 : sqrt(raiz[i]))) / 2;
        xr[i+3] = - xr[i+1];

        // Intercambio de los coeficientes para el siguiente polinomio reducido
        aux = a;
        a = b;
        b = aux;

        // Disminucion del grado del polinomio
        for (o = 0, m = 2; m <= n; m++, o++)
        {
            a[o] = a[m];
        }

        n -= 2;
    }

    if(n == 1)
    {
        xr [i] = - a [0];
    }

    return xr;
}

int main(int argc, const char * argv[])
{
    int n, i, j = 0, ni;

    float *a, *b, *db, e_a, *xr, r, s;

    printf("\nMETODO DE BAIRSTOW\n\n");

    do
    {
        printf("Ingrese el grado del polinomio\n");
        scanf("%d", &n);
    }
    while(n < 3);

    do
    {
        printf("\nIngrese el numero de iteraciones\n");
        scanf("%d", &ni);
    }
    while (ni < 1);

    do
    {
        printf("\nIngrese el error de aproximacion\n");
        scanf("%f", &e_a);
    }
    while (e_a < 0);


    a = (float*) malloc ((n + 1) * sizeof (float));

    if(a == NULL)
    {
        return 1;
    }

    b = (float*) malloc ((n + 1) * sizeof (float));

    if(b == NULL)
    {
        free(a);

        return 2;
    }

    db = (float*) calloc(n + 1, sizeof(float));

    if(db == NULL)
    {
        free(a);
        free(b);

        return 3;
    }

    xr = (float*) calloc(n * 2, sizeof(float));

    if(xr == NULL)
    {
        free(a);
        free(b);
        free(db);

        return 4;
    }

    printf("\n");

    for(i = 0; i <= n; i++)
    {
        printf("A[%d] = ", i);
        scanf("%f", a+i); // scanf("%f", &a[i]); -> &*(a+i)
    }

    printf("\nIngrese el valor inicial de r\n");
    scanf("%f", &r);

    printf("\nIngrese el valor inicial de s\n");
    scanf("%f", &s);
    
    xr = bairstow (n, a, b, db, e_a, ni, r, s);

    printf("\n");

    for(i = 0, j = 0; i < (n * 2); i+=2, j++)
    {
        printf("\nxr[%d] = (%f %+f i)", j+1, xr[i], xr[i+1]);
    }

    printf("\n\n");

    free(a);
    free(b);
    free(db);
    free(xr);

    return 0;
}
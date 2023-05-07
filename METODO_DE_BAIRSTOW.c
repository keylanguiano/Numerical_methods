#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float* bairstow (int n, float *a, float *b, float *db, float e_a, int ni,  float r, float s);

float* bairstow (int n, float *a, float *b, float *db, float e_a, int ni, float r, float s)
{
    int j;

    float R[2], vr[n], vs[n], det, dr, ds, raiz[n], ea, ear, eas, *xr;

    xr = (float*) calloc(n * 4, sizeof(float));

    // Asignación de valores iniciales para r y s
    vr[0] = r;
    vs[0] = s;

    // Ciclo que indica la continuidad del programa de acuerdo al número máximo de iteraciones
    for( j = 0, ea = 100 ; j < ni; j++)
    {
        b[n] = a[n];
        b[n-1] = a[n-1] + b[n] * vr[j];

        for(int i = n-2; i >= -1; i--)
        {
            b[i] = a[i] + b[i+1] * vr[j] + b[i+2] * vs[j];
        }

        R[1] = b[1]; // Lineal
        R[0] = b[0] - b[1] * vr[j]; // Independiente
       
        db[n]=0;
        db[n-1] = b[n];

        for(int i = n-2; i > -1; i--)
        {
            db[i] = db[i+1] * vr[j] + db[i+2] * vs[j] + b[i+1];
        }

        det = 1 / (db[2] * db[0] - db[1] * db[1]);
        
        dr = det * (db[1] * b[1] - db[2] * b[0]);
        ds = det * (db[1] * b[0] - db[0] * b[1]);
        
        vr[j+1] = vr[j] + dr;
        vs[j+1] = vs[j] + ds;

        ear = fabs (dr / vr[j+1]);
        eas = fabs (ds / vs[j+1]);

        // Selección del error de aproximación para la iteración
        if (eas > ea)
        {
            ea = eas;
        }else
        {
            ea = ear;
        }

        printf("\n%d. R[1] = %f\tR[0] = %f (%f)", j+1, R[1], R[0], ea);

        // Comprobación de la continuidad del programa por error de aproximación
        if(ea < e_a)
        {
            break;
        }

    }

    printf("\n");

    // Almacenamiento de los valores reales e imaginarios en el vector según la paridad
    for (int i = 0, j = 0; i < (n * 4); i+=4, j++)
    {
        raiz[j] = vr[j] * vr [j]+ 4 * vs[j];
        
        xr[i] = (vr[j] + (raiz[j] < 0 ? 0 : sqrt(raiz[j]))) / 2;
        xr[i+1] = raiz[j] < 0? sqrt(-raiz[j]) / 2 : 0;

        xr[i+2] = (vr[j] - (raiz[j] < 0 ? 0 : sqrt(raiz[j]))) / 2;
        xr[i+3] = - xr[i+1];
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

    xr = (float*) calloc(n * 4, sizeof(float));

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

    do
    {
        printf("\nIngrese el valor inicial de r\n");
        scanf("%f", &r);
    }
    while(r < 0);

    do
    {
        printf("\nIngrese el valor inicial de s\n");
        scanf("%f", &s);
    }
    while(s < 0);
    
    xr = bairstow (n, a, b, db, e_a, ni, r, s);

    printf("\n");

    for(i = 0, j = 0; i < (n * 4); i+=4, j+=2)
    {
        printf("\nxr[%d, %d] = (%f %+f i , %f %+f i)", j+1, j+2, xr[i], xr[i+1], xr[i+2], xr[i+3]);
    }

    printf("\n\n");

    free(a);
    free(b);
    free(db);
    free(xr);

    return 0;
}
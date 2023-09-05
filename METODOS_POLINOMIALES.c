#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void capturar(float *pa, int n);
void imprimir(float *pa, int n) ;
float evaluar2(float *pa, int n, float x);
float deflacion(float *pa, int n, float r, float *pc);
float muller(float *pa, int n, float x0, float x1, float x2, int nit, float ea);

void capturar(float *pa, int n) 
{
    int i;

    printf("\nIngrese los coeficientes\n");

    for(int i=0; i<=n; i++)
    {
        printf("A[%d] = ", i);
        scanf("%f", pa+i);
    }
}

void imprimir(float *pa, int n) 
{
    int i;
    printf("f%d(x) = %.2f", n, pa[0]);
   
    for (i = 1; i < (n + 1); i++)
        printf("% + .2fx^%d", pa[i], i);
   
    printf("\n");
}

float evaluar2(float *pa, int n, float x) 
{
    if (n)
        return evaluar2(pa + 1, n - 1, x) * x + pa[0];
    else
        return pa[0];
}

float deflacion(float *pa, int n, float r, float *pc)
{
    int i;
    pc[n - 1] = pa[n];
   
    for (i = n - 1; i > -1; i--)
        pc[i - 1] = pa[i] + pc[i] * r;
    
    return pa[0] + pc[0] * r;
}

float muller(float *pa, int n, float x0, float x1, float x2, int nit, float ea) {
    int i;
    float x[3], fx[3], h0, d0, h1, d1, a, b, c, xr1, xr2, fxr1, fxr2, raiz;
    
	x[0] = x0;
    x[1] = x1;
    x[2] = x2;
    
	for (i = 0; i < 3; i++)
        fx[i] = evaluar2(pa, n, x[i]);
    
	for (i = 0; i < nit; i++) {
	    h0 = x[(i + 1) % 3] - x[i % 3];
        d0 = (fx[(i + 1) % 3] - fx[i % 3]) / h0;
        h1 = x[(i + 2) % 3] - x[(i + 1) % 3];
        d1 = (fx[(i + 2) % 3] - fx[(i + 1) % 3]) / h1;
       
	    a = (d1 - d0) / (h1 + h0);
        b = d1 + a * h1;
        c = fx[(i + 2) % 3];
        
		raiz = b * b - 4 * a * c;
    
	    if (raiz < 0)
            return 0;
        
		xr1 = (-2 * c) / (b + sqrt(raiz));
        xr2 = (-2 * c) / (b - sqrt(raiz));
        
		x[i % 3] = x[(i + 2) % 3] + (fabs(x[(i + 2) % 3] - xr1) < fabs(x[(i + 2) % 3] - xr2) ? xr1 : xr2);
        fx[i % 3] = evaluar2(pa, n, x[i % 3]);
       
	    if (fabs((x[i % 3] - x[(i + 2) % 3]) / x[i % 3]) < ea)
            break;
    }
    return x[(i + 2) % 3];
}

void Bairstow(float *a, float *b, float *db, float error, int max_iter, int n)
{
    int i, j;
    float r=0, s=0, dr, ds, raiz, xr1r, xr1i, xr2r, xr2i, R[2], det, r0, s0;
    bool ea = false, mi = false;

    // MIENTRAS EL POLINOMIO NO SEA DE SEGUNDO GRADO SE SEGUIRA EJECUTANDO
    // LO CORRECTO SERIA CALCULAR LAS RAICES DE UN POLINOMIO DE SEGUNDO GRADO CON LA
    // FORMULA GENERAL PERO POR SIMPLICIDAD DECIDI DEJARLO ASI
    while (n >= 2)
    {
        // s0 Y r0 ME AYUDAN A CALCULAR EL ERROR DE APROXIMACION
        // EN LA PRIMER ITERACION EL VALOR ES INFINITO YA QUE AL SER EL PRIMER CICLO
        // EL ERROR DE APROXIMACION NO EXISTE YA QUE NOA HAY DATOS PREVIOS QUE NOS
        // PERMITAN HACER LA COMPARACION
        s0 = INFINITY;
        r0 = INFINITY;
        for(j=0; j<max_iter; j++)
        {
            // CALCULO DE b POR DEFLACION
            // LOS VALORES DE a PASAN A b Y SE REALIZA LA DEFLACION
            b[n] = a[n];
            b[n-1] = a[n-1]+b[n]*r;
            for(i=n-2; i>=-1; i--)
                b[i] = a[i]+b[i+1]*r+b[i+2]*s;
            // AL FINAL TENEMOS DOS ECUACION:
            // LA LINEAL
            R[1] = b[1];
            // Y LA INDEPENDIENTE
            R[0] = b[0]-b[1]*r;

            // SE HACEN LOS CALCULOS PARA OBTENER LOS VALORES DE r Y s MAS CERCANOS AL
            // ERROR DE APROXIMACION INTRODUCIDO
            db[n]=0;
            db[n-1]=b[n];
            for(i=n-2; i>-1; i--)
                db[i]=db[i+1]*r+db[i+2]*s+b[i+1];
            det = 1/(db[2]*db[0]-db[1]*db[1]);
            dr = det*(db[1]*b[1]-db[2]*b[0]);
            ds = det*(db[1]*b[0]-db[0]*b[1]);
            r+=dr;
            s+=ds;

            // EN CADA CICLO SE VALIDA EL ERROR DE APROXIMACION DE r Y s. EL PRIMERO EN
            // LLEGAR AL VALOR INTRODUCIDO ES EL QUE DETIENE EL CICLO
            if (fabs((r-r0)/r) < error || fabs((s-s0)/s) < error)
            {
                // ea ES UNA BANDERA QUE ME AYUDA A SABER SI BAIRSTOW SE DETUVO
                // PORQUE LLEGO AL ERROR DE APROXIMACION
                ea = true;
                break;
            }

            // GUARDAMOS LOS VALORES DE r Y s EN r0 Y s0 ESTO PARA CALCULAR EL ERROR DE
            // APROXIMACION EN EL SIGUIENTE CICLO
            r0 = r;
            s0 = s;

            // mi ES UNA BANDERA QUE ME AYUDA A SABER SI EL CICLO SE DETUVO PORQUE LLEGO
            // AL NUMERO MAXIMO DE ITERACIONES
            if(j+1 == max_iter)
                mi = true;
        }

        // UNA VEZ DETENIDO EL CICLO CALCULAMOS LAS RAICES USANDO LOS VALORES DE r Y s
        // EL CALCULO SE REALIZA CON LA FORMULA GENERAL
        raiz = r*r+4*s;
        
        xr1r = (r+(raiz<0?0:sqrt(raiz)))/2;
        xr1i = raiz<0?sqrt(-raiz)/2:0;
        
        xr2r = (r-(raiz<0?0:sqrt(raiz)))/2;
        xr2i = -xr1i;

        // IMPRIMIMOS LAS RAICES OBTENIDA EN LA CONSOLA
        printf("xr%d = %f%+fi\nxr%d = %f%+fi\n", n, xr1r, xr1i, n-1, xr2r, xr2i);

        // REDUCIMOS EL GRADO DEL POLINOMIO
        for (int k=n; k>1; k--)
            a[k-2] = b[k];

        // REDUCIMOS EL TAMAÑO DEL POLINOMIO
        n -= 2;
    }

    // SI n ES IGUAL A 1 SIGNIFICA QUE NUESTRO POLINOMIO ES LINEAL Y SOLO QUEDA
    // DESPEJAR EL COEFICIENTE
    if(n == 1)
    {
        a[0] *= -1;
        printf("xr1 = %f%+fi\n", a[0], b[0]);
    }

    // REVISAMOS NUESTRAS DOS BANDERAS, SI SE LLEGO AL ERROR DE APROXIMACION O
    // AL MAXIMO DE ITERACIONES O A AMBOS
    if(ea == true)
        printf("\nSE ALCANZO EL ERROR DE APROXIMACION DESEADO.\n");
    else if(mi == true)
        printf("\nSE ALCANZO EL MAXIMO DE ITERACIONES.\n");
}

int main(void) {
    float *pa, *pb, *pc, *pd, *xr, R, raiz, error, x0, x1, x2;
    int n, nb, nt, i, max_iter;

    printf("METODOS POLINOMIALES\n\n");

    do {
        printf("Ingrese el grado del polinomio (Mayor a 3)\n");
        scanf("%d", &n);
    } while (n < 3);

    do {
        printf("\nIngrese el numero de iteciones\n");
        scanf("%d", &max_iter);
    } while (max_iter < 1);

    do {
        printf("\nIngrese el error de aproximacion\n");
        scanf("%f", &error);
    } while (error < 0);

    printf("\nIngrese x0\n");
    scanf("%f", &x0);

    printf("\nIngrese x1\n");
    scanf("%f", &x1);

    printf("\nIngrese x2\n");
    scanf("%f", &x2);

    pa = (float *)malloc((n + 1) * sizeof(float));
    if (pa == NULL)
        return 1;

    xr = (float *)malloc(n * sizeof(float));
    if (xr == NULL)
        return 2;

    for(i=n, nt=1; i>2; i--)
        nt*=i;

    pd = (float *)malloc(nt*sizeof(float));
    if(pd==NULL)
        return 3;

    capturar(pa, n);

    printf("\n");
    imprimir(pa, n);

    printf("\nMuller\n");
    for (i = 0, pb = pa, pc = pd, nb = n; i < (n - 2); i++)
    {
        xr[i] = muller(pb, nb, x0, x1, x2, max_iter, error);
        
		R = deflacion(pb, nb, xr[i], pc);
        printf("Xr[%d] = %f (%f)\n", i + 1, xr[i], R);
        
		pb = pc;
        nb--;
        pc += nb;
    }

    pc = pb;
    raiz = pc[1] * pc[1] - 4 * pc[2] * pc[0];
    
    xr[n - 2] = -2 * pc[0] / (pc[1] + sqrt(raiz));
    printf("xr[%d] = %f\n", n - 1, xr[n - 2]);
    
    xr[n - 1] = -2 * pc[0] / (pc[1] - sqrt(raiz));
    printf("xr[%d] = %f\n", n, xr[n - 1]);

    // Método de Bairstow
    float *a, *b, *db;
    
    a = (float*)malloc((n+1)*sizeof(float));
    if(a==NULL)
        return 1;

    b = (float*)malloc((n+1)*sizeof(float));
    if(b==NULL)
    {
        free(a);
        return 2;
    }

    db=(float*)calloc(n+1,sizeof(float));
    if(db==NULL)
    {
        free(a);
        free(b);
        return 3;
    }

    printf("\nBairstow\n");
    Bairstow(pa, b, db, error, max_iter, n);

    free(a);
    free(b);
    free(db);
    free(pa);
    free(xr);
    free(pd);
    return 0;
}
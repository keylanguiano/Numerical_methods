#include <math.h>
#include <stdio.h>

float f(float x) { return sin(2 * exp(x / 3)); }

float biseccion(float xi, float xs, float e_a, int ni)
{
    float xm, axm, fxi, fxs, fxm, ea;
    int i;

    fxi = f(xi);
   
    if (!fxi)
        return xi;
    
    fxs = f(xs);
  
    if (!fxs)
        return xs;
   
    for (i = 0, ea = 100, xm = 100; (i < ni) && (fxi * fxs < 0); i++)
    {
        axm = xm;
        xm = (xi + xs) / 2;
        fxm = f(xm);
        
        ea = fabs((xm - axm) / xm);
        
        if (!fxm)
            return xm;
        else if (fxm * fxs < 0)
        {
            xi = xm;
            fxi = fxm;
        }
        else
        {
            xs = xm;
            fxs = fxm;
        }
        
        printf("%d. f(%f) = %f (%f)\n", i + 1, xm, fxm, ea);
        
        if (ea < e_a)
            break;
    }
    
    return xm;
}

float fposicion(float xi, float xs, float e_a, int ni)
{
    float xfp, axfp, fxi, fxs, fxfp, ea, a0, a1, d;
    int i;
   
    fxi = f(xi);
   
    if (!fxi)
        return xi;
   
    fxs = f(xs);
   
    if (!fxs)
        return xs;
   
    for (i = 0, ea = 100, xfp = 100; (i < ni) && (fxi * fxs < 0); i++)
    {
        axfp = xfp;
        d = xi - xs;
        
        a1 = fxi / d - fxs / d;
        a0 = -fxi * xs / d + fxs * xi / d;
        
        xfp = -a0 / a1;
        fxfp = f(xfp);
        
        ea = fabs((xfp - axfp) / xfp);
       
        if (!fxfp)
            return xfp;
        else if (fxfp * fxs < 0)
        {
            xi = xfp;
            fxi = fxfp;
        }
        else
        {
            xs = xfp;
            fxs = fxfp;
        }
        
        printf("%d. f(%f) = %f (%f)\n", i + 1, xfp, fxfp, ea);
        
        if (ea < e_a)
            break;
    }
    return xfp;
}

float ridder(float xi, float xs, float e_a, int ni)
{
    float xr, axr, fxi, fxs, fxr, ea, a0, a1, a2, d1, d2, d3, xr1, xr2, xm, fxm;
    int i;
   
    fxi = f(xi);
   
    if (!fxi)
        return xi;
   
    fxs = f(xs);
   
    if (!fxs)
        return xs;
   
    for (i = 0, ea = 100, xr = 100; (i < ni) && (fxi * fxs < 0); i++)
    {
        axr = xr;
        xm = (xi + xs) / 2;
        fxm = f(xm);
        
        d1 = xi * xm - xi * xs - xm * xs + xs * xs;
        d2 = xi * xm - xi * xs + xm * xs - xm * xm;
        d3 = xi * xm + xi * xs - xm * xs - xi * xi;
        
        a2 = -fxi / d3 - fxm / d2 + fxs / d1;
        a1 = (xm + xs) * fxi / d3 + (xi + xs) * fxm / d2 - (xi + xm) * fxs / d1;
        a0 = -xm * xs * fxi / d3 - xi * xs * fxm / d2 + xi * xm * fxs / d1;
        
        xr1 = (-a1 + sqrt(a1 * a1 - 4 * a2 * a0)) / (2 * a2);
        xr2 = (-a1 - sqrt(a1 * a1 - 4 * a2 * a0)) / (2 * a2);
        xr = ((xi < xr1) && (xr1 < xs)) ? xr1 : xr2;
        
        fxr = f(xr);
        
        ea = fabs((xr - axr) / xr);
        
        if (!fxr)
            return xr;
        else if (fxr * fxm < 0)
        {
            xi = xm;
            fxi = fxm;
            xs = xr;
            fxs = fxr;
        }
        else if (fxr * fxs < 0)
        {
            xi = xr;
            fxi = fxr;
        }
        else
        {
            xs = xr;
            fxs = fxr;
        }
        
        printf("%d. f(%f) = %f (%f)\n", i + 1, xr, fxr, ea);
        
        if (ea < e_a)
            break;
    }
    return xr;
}

int main(void)
{
    float xi, xs, xrb, xrfp, xrr, e_a;
    int ni;
    
    printf("METODOS CERRADOS\n\n");

    printf("Ingrese el valor inferior\n");
    scanf("%f", &xi);
    
    printf("Ingrese el valor superior\n");
    scanf("%f", &xs);
    
    if (xi > xs)
    {
        xi *= xs;
        xs = xi / xs;
        xi /= xs;
    }
    
    do
    {
        printf("Ingrese el numero de iteraciones\n");
        scanf("%d", &ni);
    } while (ni < 1);
   
    do
    {
        printf("Ingrese el error de aproximacion\n");
        scanf("%f", &e_a);
    } while (e_a < 0);
    
    printf("Metodo de Biseccion:\n");
    xrb = biseccion(xi, xs, e_a, ni);
    
    printf("Metodo de Falsa Posicion:\n");
    xrfp = fposicion(xi, xs, e_a, ni);
    
    printf("Metodo de Ridder:\n");
    xrr = ridder(xi, xs, e_a, ni);
    
    return 0;
}
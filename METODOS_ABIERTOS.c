#include <stdio.h>
#include <math.h>

#define Dx 1e-5

float f(float x)
{
    return 2 * pow(x, 3) - 11.7 * pow(x, 2) + 17.7 * x - 5;
}

float g(float x)
{
    return f(x) + x;
}

float df(float x)
{
    return 6 * pow(x, 2) - 23.4 * x + 17.7;
}

float pf(float xi, float e_a, int ni)
{
  float xpf, axpf, ea, fxpf;
  int i;
  xpf = xi;
  fxpf = f(xpf);
  if(!fxpf)
    return xpf;
  for(i=0, ea=100; i<ni; i++)
  {
    axpf = xpf;
    xpf = g(xpf);
    ea = fabs((xpf-axpf)/xpf);
    fxpf = f(xpf);
    if(!fxpf)
      return xpf;
    printf("%d. f(%f) = %f (%f)\n", i+1, xpf, fxpf, ea);
    if(ea<e_a)
      break;
  }
  return xpf;
}

float nr(float xi, float e_a, int ni)
{
  float xnr, axnr, ea, fxnr, dfxnr;
  int i;
  xnr = xi;
  fxnr = f(xnr);
  if(!fxnr)
    return xnr;
  for(i=0, ea=100; i<ni; i++)
  {
    axnr = xnr;
    dfxnr = df(xnr);
    xnr = xnr-fxnr/dfxnr;
    ea = fabs((xnr-axnr)/xnr);
    fxnr = f(xnr);
    if(!fxnr)
      return xnr;
    printf("%d. f(%f) = %f (%f)\n", i+1, xnr, fxnr, ea);
    if(ea<e_a)
      break;
  }
  return xnr;
}

float secante(float xi, float e_a, int ni)
{
  float xs, axs, ea, fxs, dfxs;
  int i;
  xs = xi;
  fxs = f(xs);
  if(!fxs)
    return xs;
  for(i=0, ea=100; i<ni; i++)
  {
    axs = xs;
    dfxs = (f(xs+Dx)-fxs)/Dx;
    xs = xs-fxs/dfxs;
    ea = fabs((xs-axs)/xs);
    fxs = f(xs);
    if(!fxs)
      return xs;
    printf("%d. f(%f) = %f (%f)\n", i+1, xs, fxs, ea);
    if(ea<e_a)
      break;
  }
  return xs;
}

int main(void) {
  float xi, xrpf, xrnr, xrs, e_a;
  int ni;
  printf("Ingrese el valor inicial: ");
  scanf("%f", &xi);
  do{
    printf("Ingrese el numero de iteraciones:");
    scanf("%d", &ni);
  }while(ni<1);
  do{
    printf("Ingrese el error de aproximacion: ");
    scanf("%f", &e_a);
  }while(e_a<0);
  printf("Metodo de punto fijo:\n");
  xrpf = pf(xi, e_a, ni);
  printf("Metodo de Newton-Raphson:\n");
  xrnr = nr(xi, e_a, ni);
  printf("Metodo de la Secante:\n");
  xrs = secante(xi, e_a, ni);
  return 0;
}
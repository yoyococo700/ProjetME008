#include <stdio.h>
#include <stdlib.h>
#define n 100

void factoriser_tridiago(float d[n], float c[n], float a[n], float l[n],
                         float u[n], float v[n]) {
  int i;
  c[n - 1] = 0, v[n - 1] = 0, l[0] = 0;
  u[0] = d[0];
  for (i = 0; i < n; i++) {
    v[i] = c[i];
    l[i + 1] = a[i + 1] / u[i];
    u[i + 1] = d[i + 1] - l[i + 1] * v[i];
  }
}

void afficher_vect(float x[n]) {
  int i;
  printf("--------------------------------------\n");
  for (i = 0; i < n; i++)
    printf("%f\n", x[i]);
  printf("--------------------------------------\n");
}

void CopyTab(float a[n], float b[n]) {
  for (int i = 0; i < n; i++)
    a[i] = b[i];
}

void init(float h, float a[n], float c[n], float d[n]) {
  for (int i = 0; i < n; i++) {
    a[i] = h;
    c[i] = h;
    d[i] = 1 - 2 * h;
  }
  d[0] = 1 - h;
  d[n - 1] = 1 - h;
}

void initU(float Un[n]) {
  for (int i = 0; i < n; i++) {
    if (i < 50)
      Un[i] = 1;
    else
      Un[i] = 10;
  }
  Un[49] = 5.5;
}

void updateUn(float a[n], float c[n], float d[n], float Un[n]) {
  float Un1[n];
  Un1[0] = d[0] * Un[0] + c[0] * Un[1];
  Un1[n - 1] = d[n - 1] * Un[n - 1] + a[n - 1] * Un[n - 2];
  for (int i = 1; i < n - 1; i++)
    Un1[i] = a[i] * Un[i - 1] + d[i] * Un[i] + c[i] * Un[i + 1];

  CopyTab(Un, Un1);
}

void Un_instant_t(float Unt[n], float instant, float dt, float dx) {
  float a[n], c[n], d[n];
  float l[n], u[n], v[n];
  float h = dt / (dx * dx);
  init(h, a, c, d);
  initU(Unt);
  factoriser_tridiago(d, c, a, l, u, v);
  for (int i = 0; i * dt < instant; i++) {
    updateUn(a, c, d, Unt);
    }
}

int main(int argc, char *argv[]) {

  float dt, dx;
  scanf("%f", &dt);
  scanf("%f", &dx);


  float  Un0[n],Un1[n], Un2[n], Un3[n], Un4[n], Un20[n];


//Question 2
  // float a[n], c[n], d[n];
  // float l[n], u[n], v[n];
  // float h = dt / (dx * dx);
  // init(h, a, c, d);
  // initU(Un0);
  // factoriser_tridiago(d, c, a, l, u, v);
  // for (int i = 0; i*dt < 200; i++)
  // {
  //     updateUn(a,c,d,Un0);
  //     afficher_vect(Un0);
  // }
  //Fin de la question 2


  //Question 3
  Un_instant_t (Un1, 1, dt, dx);
  Un_instant_t (Un2, 2, dt, dx);
  Un_instant_t (Un3, 3, dt, dx);
  Un_instant_t (Un4, 4, dt, dx);
  Un_instant_t (Un20, 20, dt, dx);

  FILE *out1 = fopen("question3.txt","wt");
  for(int i = 0; i<n; i++) {
    fprintf(out1,"%f\t%f\t%f\t%f\t%f\t%f\n", i*dx, Un1[i],Un2[i],Un3[i],Un4[i],Un20[i]);
  }
  fclose(out1);
  system("gnuplot question3.p");
 //Fin de la question3 
  
  
  return 0;
}
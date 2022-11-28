#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define n 100
#define TAILLE_INTERVALLE 20


void diff_vect(float U1[n],float U2[n], float res[n]){
  for (int i = 0; i < n; i++)
    res[i]=U2[i]-U1[i];
}

float norme1(float U[n]){
  float s = 0;
  for (int i = 0; i < s; i++)
    s+=U[i];
  return s;  
}

float interpolation_lineaire(float t, float dt2){
  float proportionalite = TAILLE_INTERVALLE/(n*1.0);

}


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

void question2(float dt,float dx){
  float  Un0[n];

  //Question 2
  float a[n], c[n], d[n];
  float l[n], u[n], v[n];
  float h = dt / (dx * dx);
  init(h, a, c, d);
  initU(Un0);
  factoriser_tridiago(d, c, a, l, u, v);
  for (int i = 0; i*dt < 200; i++)
  {
      updateUn(a,c,d,Un0);
      afficher_vect(Un0);
  }
  //Fin de la question 2

}

void enregistrer_vect(float tmax, float dt, float dx, char titre[10] ){

  float  Un0[n];
  FILE *out1 = fopen(titre,"wt");
  float a[n], c[n], d[n];
  float l[n], u[n], v[n];
  float h = dt / (dx * dx);
  init(h, a, c, d);
  initU(Un0);
  factoriser_tridiago(d, c, a, l, u, v);
  for (int i = 0; i*dt < tmax; i++)
  {
      updateUn(a,c,d,Un0);
      fprintf(out1,"%f\t",i*dt);
      for (int j = 0; j < n; j++)
      {
        fprintf(out1,"%f\t",Un0[j]);
      }
      fprintf(out1,"\n");
  }
  fclose(out1);
}

void question3(float dt,float dx){
  float  Un1[n], Un2[n], Un3[n], Un4[n], Un20[n];
  
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

}

void question3_2(float dt,float dx){
  int m = 5;
  float liste_points[5] = {-8,-4,0,4,8};

  float  Unt[n];
  float a[n], c[n], d[n];
  float l[n], u[n], v[n];
  float h = dt / (dx * dx);
  float proportionalite = TAILLE_INTERVALLE/(n*1.0);
  
  init(h, a, c, d);
  initU(Unt);
  factoriser_tridiago(d, c, a, l, u, v);

  FILE *out1 = fopen("question3_2.txt","wt");
  
 
  for (int i = 0; i * dt < 200; i++) {
    fprintf( out1, "%f\t", i*dt); //Impression du temps au debut de la ligne 
    for (int j = 0; j < m; j++)
    {
      int k = ((int)(liste_points[j]/proportionalite))+(n/2);
      fprintf(out1, "%f\t",Unt[k]);
    }
    fprintf(out1,"\n");
    
    updateUn(a, c, d, Unt);
  }
  fclose(out1);
  system("gnuplot question3_2.p");
}

float compare_dt(float tmax, float dt1, float dt2,float dx){
  //dt1 doit toujours etre inferieur a dt2 pour la suite
  if (dt1>dt2){
    float temp = dt2;
    dt2= dt1;
    dt1 = temp;
  }

  float h = dt2/dt1;
  float t1,t2;  
  float U1[n],U2[n];


  enregistrer_vect(tmax,dt1,dx,"1.dat");
  enregistrer_vect(tmax,dt2,dx,"2.dat");
  
  FILE* in = fopen("2.dat","r");
  FILE* out = fopen("inter2.dat","wt");


  //On va maintenant interpoller 2 sur 1
  for (int i = 0; i*dt2 < tmax; i++)
  { 
    fscanf(in,"%f",&t1);
    for (int j = 0; j<n; j++)
      fscanf(in,"%f",&U1[j]);
    fscanf(in,"%f",&t2);
    for (int j = 0; j<n; j++)
      fscanf(in,"%f",&U2[j]);
    printf("%f\t%f\n",t1,t2);
    for (int j = t1; j < h; j++)
    {
      
    }
    
  }
  
  fclose(in);

}


int main(int argc, char *argv[]) {

  //question2(0.001,1);
  // question3(0.001,1);
  // question3_2(0.001,1);
  compare_dt(200,0.001,0.001,1);
  
  
  return 0;
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define n 100
#define TAILLE_INTERVALLE 20


typedef struct sim_param{
  int N;
  float taille_intervalle;
  float tmax;
  float dt;
  float dx;
}sim_param;


// typedef struct Unt{
//   sim_param param;
//   float t;
//   float *data;
// }Unt;

typedef struct saved_V{
  char name[10];
  sim_param param;
  float** data;
}saved_V;


void affiche_sim_param(sim_param p){
  printf("---------------Parametres de simulation-----------------\n");
  printf("Echantillonage spatial N = %d\n",p.N);
  printf("Duree de la simulation = %f\n",p.tmax);
  printf("Largeur de la simulation = %f\n",p.taille_intervalle);
  printf("dt = %f\n",p.dt);
  printf("dx= %f\n",p.dx);
  printf("nb_lignes = %d\n",(int)(p.tmax/p.dt)+1);
  printf("nb_colonne = %d\n", p.N);
  printf("---------------Fin des parametres de simulation-----------------\n");
}

sim_param standard_param(float dt,float dx){
  sim_param param;
  param.dt = dt;
  param.dx = dx;
  param.N = n;
  param.taille_intervalle = TAILLE_INTERVALLE;
  param.tmax = 200;
  return param;
}

int compare_sim_param(sim_param param1,sim_param param2){
  if (param1.dt!=param2.dt)
    return 0;
  if (param1.dx!=param2.dx) 
    return 0;
  if (param1.N!=param2.N)
    return 0;
  if (param1.tmax!=param2.tmax)
    return 0;
  if (param1.taille_intervalle!=param2.taille_intervalle)
    return 0;
  return 1;
}

void affiche_saved_V(saved_V* V){
  printf("Debut de affiche save \n");
  for (int i = 0; i*V->param.dt < V->param.tmax; i++)
  { 
    printf("%d\t",i);
    printf("%f:\t",V->data[i][V->param.N]);
    for (int j = 0; j < V->param.N; j++)
    {
      printf("%f\t",V->data[i][j]);
    }
    printf("\n");
  }  
}


saved_V* create_empty_savedV(sim_param param){

  int lenght = param.N;
  int nb_lines = (int)(param.tmax/param.dt)+1;
  saved_V* res = (saved_V*)malloc(sizeof(saved_V));

  if (res==NULL)
  {
    printf("erreur dans la creation du vecteur\n");
    exit;
  }
  
  res->param = param;
  res->data = (float**)malloc(nb_lines*sizeof(float*));

  if (res->data==NULL)
  {
    printf("erreur dans la creation du tableau\n");
    printf("L'erreur provient peut etre du nombre de ligne = %d\n",nb_lines);
    exit;
  }
 
  for (int i = 0; i < nb_lines; i++){
    res->data[i] = (float*)malloc((lenght+1)*sizeof(float));
    if (res->data[i]==NULL)
    {
      printf("erreur dans la creation du vecteur dans le tableau\n");
      exit;
    }
    
  }

  for (int i = 0; i < nb_lines; i++)
    for (int j = 0; j < lenght+1; j++)
      res->data[i][j]=0;
  
  
  //printf("%f\n",res->data[5][1]);

  printf("Initialisation du saved_V reussi\n");
  return res;
}

void free_savedV(saved_V* V){
  int lenght = V->param.N;
  int nb_lines = (int)(V->param.tmax/V->param.dt);
  for (int i = 0; i < lenght+1; i++)
  {
    free(V->data[i]);
  }
  free(V->data);
  free(V);
  printf("Destruction du vecteur reussi\n");
}

//Unt[n] = t
void save_Unt_in_saved_V(saved_V* V,float Unt[],float t,sim_param param){
  if (compare_sim_param(param,V->param))
  { 
    int k = (int)(t/param.dt);
    //printf("%d\n",k);              // k correspond a la ligne d'insertion 
    if (V->data[k]==NULL)
    {
      printf("attention, depassement d'indice k=%d\nSortie de la fonction save_unt\n",k);
      return;
    }
    
    for (int i = 0; i < param.N; i++){
      V->data[k][i]=Unt[i];
    }
    
    V->data[k][param.N] = t;
  }
  else{
    printf("Le vecteur entre ne correspond pas au format desiree\n");
    printf("Le vecteur n'as pas pu etre enregistre\n");
  }
}

void diff_vect(float U1[],float U2[], float res[],sim_param p){
  for (int i = 0; i < p.N; i++)
    res[i]=U2[i]-U1[i];
}

float norme1(float U[],sim_param p){
  float s = 0;
  for (int i = 0; i < p.N; i++)
    s+=fabs(U[i]);
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

saved_V* enregistrer_vect_V(sim_param param){
  float dt = param.dt;
  float dx =param.dx;
  float tmax = param.tmax;
  float  Un0[n];

  
  float a[n], c[n], d[n];
  float l[n], u[n], v[n];
  float h = dt / (dx * dx);

  saved_V *res = create_empty_savedV(param);
  
  //affiche_sim_param(param);

  init(h, a, c, d);
  initU(Un0);

  for (int i = 0; i*dt < tmax; i++)
  {
      updateUn(a,c,d,Un0);
      //printf("debug enregistrer_vect_V %d\n",i);
      save_Unt_in_saved_V(res,Un0,i*dt,param);
      
  }
  return res;
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

float compare_dt(float dt1, float dt2,float dx){
  //dt1 doit toujours etre inferieur a dt2 pour la suite
  //dt1<dt2
  if (dt1>dt2){
    float temp = dt2;
    dt2= dt1;
    dt1 = temp;
  }
  sim_param p1 = standard_param(dt1,dx);
  sim_param p2 = standard_param(dt2,dx);

  int k = (int)(dt2/dt1);
 
  float* res = (float*)malloc((p1.N+1)*sizeof(float));
  float s=0; 



  saved_V* V1 = enregistrer_vect_V(p1);
  saved_V* V2 = enregistrer_vect_V(p2);

  for (int i = 0; i*dt2 < p2.taille_intervalle ; i++)
  {
    diff_vect(V1->data[k*i],V2->data[i],res,p2);
    s+=norme1(res,p1);
  }   
  
  free_savedV(V1);
  free_savedV(V2);
  free(res);
  return s;

}


void question4(float dt,float dx){

  FILE* out = fopen("question4.dat","wt");
  if (out == NULL)
  {
    printf("erreur dans la creation du fichier de sortie\n");
    return;
  }
  
  ;

  for (int i = 4; i < 100; i+=2){
    float s = compare_dt(dt,i*dt,dx);
    fprintf(out,"%d\t%f\n",i,s);
    
  }
  fclose(out);
  system("gnuplot question4.p");
}


int main(int argc, char *argv[]) {


  //printf("%f\n",compare_dt(0.001,2,1));
  //question2(0.001,1);
  //question3(0.001,1);
  // question3_2(0.001,1);
  question4(0.001,1);
  
  return 0;
}
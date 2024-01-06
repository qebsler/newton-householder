#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#define EPSILON 0.001
#define MAX_ITER 100
using namespace std;

// [a_0, a_1, ..., a_n] -> P = a_0 + a_1*X + ... + a_n*X^n
typedef vector<complex<double>> polynome;

// Génère un complexe au hasard dans la boule (norme 1) de rayon sup
double rand(const double &sup) {
  long prec = 10000;
  return -sup + 2 * sup * (random() % prec) / prec;
}

// Calcul la dérivée d'un polynome
polynome derivee(const polynome &v){
    int l = v.size();
    polynome d = {};
    for (double i=1; i<l; i++){
        d.push_back(v[i]*i);
    }
    return d;
}

// Retourne le degré d'un polynome
int deg(const polynome &p) {
  return p.size() - 1;
}

// Calcul P(z) avec la méthode de Horner
complex<double> eval(const polynome &p, const complex<double> &z) {
  int d = deg(p);
  complex<double> coeff = p[d];

  for (int k = d - 1; k >= 0; k--) coeff = coeff * z + p[k];

  return coeff;
}

// Retourne la distance à la racine la plus proche
double dist_root(const vector<complex<double>> &roots, const complex<double> &z) {
  double m = abs(roots[0] - z);
  for (complex<double> root: roots) {
    double d = abs(root - z);
    if (d < m) m = d;
  }

  return m;
}

// Retourne l'indice de la racine la plus proche dans le vecteur roots
int id_nearest_root(const vector<complex<double>> &roots, const complex<double> &z) {
  double m = abs(roots[0] - z);
  int i = 0;

  for (int k = 1; k < roots.size(); k++) {
    double d = abs(roots[k] - z);
    if (d < m) {
      i = k;
      m = d;
    }
  }

  return i;
}

// Methode de Newton pour trouver les racines du polynome
complex<double> newton (const vector<complex<double>> &rac, const complex<double> &z0, const polynome &p) {
    complex<double> z = z0;
    polynome d = derivee(p);
    int i = 0;
    while ((dist_root(rac,z) > EPSILON) && (i < MAX_ITER)) {
        z = z - eval(p, z) / eval(d, z);
        i++;
    }
    return z;
}

// Methode de Householder pour trouver les racines du polynome
complex<double> householder (const vector<complex<double>> &roots, const complex<double> &z0, const polynome &p) {
  complex<double> z = z0;
  complex<double> h, v0, v1;
  polynome p1 = derivee(p);
  polynome p2 = derivee(p1);
  int i = 0;

  while ((dist_root(roots, z) > EPSILON) && (i < MAX_ITER)) {
    v0 = eval(p, z);
    v1 = eval(p1, z);
    h = (eval(p, z) * eval (p2, z)) / (2. * v1 * v1);
    z = z - (1. + h) * v0 / v1;
  }

  return z;
}

// Methode de Newton pour générer les données graphiques
vector<double> newton_verbose (const vector<complex<double>> &roots, const complex<double> &z0, const vector<polynome> &p) {
  complex<double> z = z0, z_last;
  int id_root, i;
  vector<double> result = {real(z0), imag(z0), -1, MAX_ITER};

  for (i = 0; i < MAX_ITER; i++) {
    z_last = z;
    z = z - eval(p[0], z) / eval(p[1], z);
    if (abs(z - z_last) < EPSILON) break;
  }

  if (i < MAX_ITER) {
    id_root = id_nearest_root(roots, z);
    result[2] = (double) id_root;
    result[3] = (double) i;
  }

  return result;
}

// Méthode de Householder pour générer les données graphiques
vector<double> householder_verbose (const vector<complex<double>> &roots, const complex<double> &z0, const vector<polynome> &p) {
  complex<double> z = z0, z_last;
  int id_root, i;
  vector<double> result = {real(z), imag(z), -1, MAX_ITER};
  complex<double> h, v0, v1;

  for (i = 0; i < MAX_ITER; i++) {
    z_last = z;
    v0 = eval(p[0], z);
    v1 = eval(p[1], z);
    h = (v0 * eval(p[2], z)) / (2. * v1 * v1);
    z = z - (1. + h) * v0 / v1;
    if (abs(z - z_last) < EPSILON) break;
  }

  if (i < MAX_ITER) {
    id_root = id_nearest_root(roots, z);
    result[2] = (double) id_root;
    result[3] = (double) i;
  }

  return result;
}

// Affiche un polynome a coeff complexes
void print_pol (const polynome &p) {
    for (int i=0; i < p.size(); i++)
        printf("(%.3f, %.3f) ", real(p[i]), imag(p[i]));
    printf("\n");
}

// p * q -> retourne le produit de p et q
polynome operator* (const polynome &p, const polynome &q) {
  polynome result = {};
  int deg_p = deg(p), deg_q = deg(q);

  for (int n = 0; n <= deg_p + deg_q; n++) {
    complex<double> a_n (0., 0.);
    for (int k = 0; k <= n; k++) {
      a_n += p[k] * q[n-k];
    }

    result.push_back(a_n);
  }

  return result;
}

// Division par des polynomes de degré 1 par la méthode de Horner
// Divise p par (X - z)
polynome rm_root(const polynome &p, const complex<double> &z) {
  int d = deg(p);
  polynome q (p.size() - 1);

  q[d - 1] = p[d];
  for (int i = d - 1; i > 0; i--) {
    q[i - 1] = q[i] * z + p[i];
  }

  return q;
}

// Trouve une racine de p par Newton
complex<double> find_root_with_z0(const complex<double> &z0, const polynome &p) {
  complex<double> z = z0, last_z;
  polynome d = derivee(p);
  double dist;
  int i = 1;

  do {
    last_z = z;
    z = z - eval(p, z) / eval(d, z);
    dist = abs(z - last_z);
    i++;
  } while (i < MAX_ITER && dist > EPSILON);

  if (i == MAX_ITER) throw invalid_argument("converge pas");
  if (!isfinite(real(z)) || !isfinite(imag(z))) throw invalid_argument("nan or inf");
  else return z;
}

// Calcul une majoration du rayon de convergence
double conv_radius(const polynome &p) {
  double m = 0;
  for (int i = 0; i < p.size(); i++) {
    if (abs(p[i]) > m) m = abs(p[i]);
  }

  return m + 1;
}

// Calcul une des racines de p
complex<double> find_root(const polynome &p) {
  complex<double> z = 0.;
  complex<double> root;
  double sup = conv_radius(p);
  bool ok = false;

  while (!ok) {
    try {
      root = find_root_with_z0(z, p);
      ok = true;
    } catch (invalid_argument) {
      complex<double> z1 (rand(sup), rand(sup));
      z = z1;
    }
  }

  return root;
}

// Calcul de toutes les racines de p
vector<complex<double>> find_roots(const polynome &p) {
  vector<complex<double>> roots = {};
  polynome q = p;
  while (deg(q) > 0) {
    complex<double> r0 = find_root(q);

    roots.push_back(r0);
    q = rm_root(q, r0);
  }

  return roots;
}

// radius en norme 1
// Génère les données graphiques
template<typename method>
void explore(const polynome &p, const vector<complex<double>> &roots, const method &m, const double &radius, const int &prec, const char* filename) {
  double step = 2 * radius / (prec);
  complex<double> step1 (step, 0), step2 (0, step);

  polynome d = derivee(p);
  vector<polynome> P = {p, d, derivee(d)};
  vector<double> reader;

  char s[20];
  sprintf(s, "data/%s.dat", filename);
  FILE* file = fopen(s, "w");

  for (double r = -radius; r <= radius; r += step) {
    for (double i = -radius; i <= radius; i += step) {
      complex<double> z (r, i);
      reader = m(roots, z, P);
      fprintf(file, "%.4f %.4f %.1f %.1f\n", reader[0], reader[1], reader[2], reader[3]);
    }
  }

  fclose(file);
}

// Lit un polynome depuis un fichier source
polynome read_pol(const char filename[]) {
  FILE* file = fopen(filename, "r");
  double real, imag;
  polynome p = {};

  if (file == NULL) throw invalid_argument("le fichier existe pas.");

  while (fscanf(file, "%lf %lf", &real, &imag) != EOF) {
    complex<double> z (real, imag);
    p.push_back(z);
  }

  fclose(file);
  return p;
}

int main(int argc, char* argv[]) {
  srand(time(NULL));
  polynome p = read_pol((char*) "polynome");
  polynome roots = find_roots(p);
  double radius;
  int resolution;

  if (argc < 3) throw invalid_argument("arguments missing");
  else {
    radius = atof(argv[1]);
    resolution = atoi(argv[2]);
  }

  print_pol(roots);

  explore(p, roots, newton_verbose, radius, resolution, (char*) "newton");
  explore(p, roots, householder_verbose, radius, resolution, (char*) "householder");

  return 0;
}


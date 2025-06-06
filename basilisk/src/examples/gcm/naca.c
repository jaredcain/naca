// CC99='mpicc -std=c99' qcc -autolink -Wall -O2 -D_MPI=1 naca.c -o naca -lm -lfb_tiny
// mpirun -np 4 ./naca <<NACA>> <<AoA>> <<Re>> <<LEVEL>>
// tail -f 0012_2.log
// http://basilisk.fr/three.js/editor/index.html?ws://COMECH-2065.ds.sc.edu:7100

#include "ibm-gcm.h"
#include "my-centered.h"
#include "ibm-gcm-events.h"
#include "view.h"

double mm=0.0, pp=0.0, tt=0.12;
double aoa = 2.0 * M_PI / 180.0;
double Re = 10000;
int maxlevel = 13;
int minlevel = 7;
const char* nacaset="0012";
const double t_end = 30;
const double framerate = 1.0/25.0;
const double chord = 1.0;
const double uref = 1.0;
const coord ci = {0.5, 0.0};
const coord cr = {0.25*chord, 0.0};

face vector muv[];
scalar omega[],m[];

void nacaset_parse(const char* nacaset, double* mm, double* pp, double* tt) {
    *mm = (nacaset[0] - '0') * 0.01;
    *pp = (nacaset[1] - '0') * 0.1;
    *tt = ((nacaset[2] - '0') * 10 + (nacaset[3] - '0')) * 0.01;
}

void naca (scalar c, face vector f, double aoa, vertex scalar phii = {0}) {
  vertex scalar phi = automatic (phii);
  foreach_vertex() {
    double xx = cr.x + (x - ci.x)*cos(aoa) - (y - ci.y)*sin(aoa);
    double yy = cr.y + (x - ci.x)*sin(aoa) + (y - ci.y)*cos(aoa);
    double xc = clamp(xx / chord, 0., 1.);
    double yc = yy / chord;
    if (xc >= 0. && xc <= 1.) {
      double yt = 5.0 * tt * (
        0.2969 * sqrt(xc) -
        0.1260 * xc -
        0.3516 * sq(xc) +
        0.2843 * cube(xc) -
        0.1015 * pow(xc, 4.0));
      double yc_camber, dyc_dx;
      if (xc < pp && pp > 1e-6) {
        yc_camber = mm / sq(pp) * (2. * pp * xc - sq(xc));
        dyc_dx = 2. * mm / sq(pp) * (pp - xc);
      } else if (pp < 1.0 - 1e-6) {
        yc_camber = mm / sq(1. - pp) * (1. - 2. * pp + 2. * pp * xc - sq(xc));
        dyc_dx = 2. * mm / sq(1. - pp) * (pp - xc);
      } else {
        yc_camber = 0.;
        dyc_dx = 0.;
      }
      phi[] = sq(yc - yc_camber) - sq(yt);
      if (xc > 0.98) {
        double d = xc - 0.98;
        phi[] += 10. * d * d;
      }
    }
    else {
      phi[] = 1.;
    }
  }
  boundary ({phi});
  fractions (phi, c, f);
}

u.n[left] = dirichlet (1);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (0)

event properties (i++) {  
  foreach_face()
    muv.x[] = fm.x[]/Re;
 }
 
event init (t = 0) {
  scalar ibmr[];
  astats ss;
  int ic = 0;
  do {
    ic++;
    naca (ibm, ibmf, aoa);
	foreach()
      ibmr[] = ibm[];
    ss = adapt_wavelet ({ibmr}, (double[]) {1.e-30}, maxlevel, minlevel);
  } while ((ss.nf || ss.nc) && ic < 100);
  naca (ibm, ibmf, aoa);
  foreach()
    u.x[] = ibm[];
  boundary((scalar *){u});
}

event logfile (i++; t <= t_end) {
  coord Fp, Fmu;
  ibm_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq((uref))*(chord));
  double CL = (Fp.y + Fmu.y)/(0.5*sq((uref))*(chord));
  vorticity (u , omega);
  stats pr  = statsf(p);
  stats ux = statsf(u.x);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", i, t, ux.max, CL, CD, pr.max, pr.min);
  fflush(stderr);
}

event movies (t += framerate, t <= t_end) {
  view(fov = 1, tx = -0.05, width = 1920, height = 1080);
  draw_vof ("ibm", "ibmf", filled = -1);
  stats pr  = statsf(p);
  squares(color = "p", max = pr.max, min = pr.min); 
  cells();
  char video[50];
  snprintf(video, sizeof(video), "%s_%.0f_%.0f_%d.mp4", nacaset, aoa * 180.0 / M_PI, Re, maxlevel);
  save(video);
}

/* event movies (t += framerate, t <= t_end) {
  view(quat = {0.0, 0.0, 0.0, 1.0}, fov = 30, near = 0.01, far = 1000, tx = -0.35, ty = 0.0, tz = -2.25, width = 1080, height = 1080);
  box();
  cells();
  draw_vof ("ibm", "ibmf", filled = -1);
  stats pr  = statsf(p);
  squares(color = "p", max = pr.max, min = pr.min, cbar = true, border = true, pos = {0.4, 0.3}, label = "p", mid = true, format = " %0.1e", levels = 100, size = 20 , lw=1, fsize = 75); 
  char video[50];
  snprintf(video, sizeof(video), "%s_%.0f_%.0f_%d.mp4", nacaset, aoa * 180.0 / M_PI, Re, maxlevel);
  save(video);
} */

event adapt (i++) {
  scalar ibmsf[];
  foreach()
    ibmsf[] = vertex_average(point, ibm);
  adapt_wavelet ({ibmsf,u}, (double[]){1.e-15,3e-3,3e-3}, maxlevel, minlevel);
}

int main(int argc, char *argv[]) {
  if (argc > 1) {
    nacaset = argv[1];
    if (argc > 2) aoa = atof(argv[2]) * M_PI / 180.0;
    if (argc > 3) Re = atof(argv[3]);
    if (argc > 4) maxlevel = atoi(argv[4]);
  }
  nacaset_parse(nacaset, &mm, &pp, &tt);
  char log_filename[50];
  snprintf(log_filename, sizeof(log_filename), "%s_%.0f_%.0f_%d.log", nacaset, aoa * 180.0 / M_PI, Re, maxlevel);
  freopen(log_filename, "w", stderr);
  L0 = 16.;
  N = 1 << (maxlevel-6);
  origin (-L0/8, -L0/2.);
  TOLERANCE = 1.e-6;
  mu = muv;
  run(); 
}

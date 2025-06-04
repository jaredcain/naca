// CC99='mpicc -std=c99' qcc -autolink -Wall -O2 -D_MPI=1 naca.c -o naca -lm -lfb_tiny
// mpirun -np 4 ./naca <<NACA>> <<AoA>> <<Re>> <<LEVEL>>
// tail -f 0012_2.log
// http://basilisk.fr/three.js/editor/index.html?ws://COMECH-2065.ds.sc.edu:7100

#include "ibm-gcm.h"
#include "my-centered.h"
#include "ibm-gcm-events.h"
#include "view.h"

double mm=0., pp=0., tt=0.12; // camber,location,thickness
double Reynolds = 10000;
double aoa = 0. * M_PI / 180.0;
const double t_end = 30;

const char* nacaset="0012";
int maxlevel = 11;
face vector muv[];

#define chord   (1.)
#define uref    (1.)

coord ci = {0.5, 0};
coord cr = {0.25*chord, 0.};

void nacaset_f(const char* nacaset, double* mm, double* pp, double* tt)
{
    *mm = (nacaset[0] - '0') * 0.01;
    *pp = (nacaset[1] - '0') * 0.1;
    *tt = ((nacaset[2] - '0') * 10 + (nacaset[3] - '0')) * 0.01;
}

int main(int argc, char *argv[]) {
  if (argc > 1) {
    nacaset = argv[1];
    if (argc > 2) aoa = atof(argv[2]) * M_PI / 180.0;
    if (argc > 3) Reynolds = atof(argv[3]);
    if (argc > 4) maxlevel = atoi(argv[4]);
  }
  
  nacaset_f(nacaset, &mm, &pp, &tt);
  
  char log_filename[50];
  snprintf(log_filename, sizeof(log_filename), "%s_%.0f_%.0f_%d.log", nacaset, aoa * 180.0 / M_PI, Reynolds, maxlevel);
  freopen(log_filename, "w", stderr);
  
  L0 = 16.;
  N = 1 << (maxlevel-4);
  origin (-L0/8, -L0/2.);
  TOLERANCE = 1.e-5 [*];
  mu = muv;
  
  run(); 
}

event properties (i++) {
  
  foreach_face()
    muv.x[] = fm.x[]/Reynolds;
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

void naca (scalar c, face vector f, double aoa, vertex scalar phii = {0})
{
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
        0.1015 * pow(xc, 4.0)
      );

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

event init (t = 0) {
  scalar ibm1[];
  astats ss;
  int ic = 0;
  do {
    ic++;
    naca (ibm, ibmf, aoa);
	foreach()
      ibm1[] = ibm[];
    ss = adapt_wavelet ({ibm1}, (double[]) {1.e-30},
                        maxlevel, minlevel = maxlevel-4);
  } while ((ss.nf || ss.nc) && ic < 100);

  naca (ibm, ibmf, aoa);
  foreach()
    u.x[] = ibm[];
  boundary((scalar *){u});
}

scalar pid[];
scalar omega[];

event logfile (i++; t <= t_end) {
	
  coord Fp, Fmu;
  ibm_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq((uref))*(chord));
  double CL = (Fp.y + Fmu.y)/(0.5*sq((uref))*(chord));

  vorticity (u , omega);
  foreach()
    pid[] = pid();
  
  stats om  = statsf(omega);
  fprintf (stderr, "%d %g %g %g %g\n", i, t, om.max, CL, CD);
  fflush(stderr);
}

event snapshot (t += 5; t <= 30) {
  char fname[80];
  sprintf(fname, "snapshot_%.0f.png", t);

  scalar omega_clean[];
  foreach()
    omega_clean[] = (ibm[] >= 1.) ? omega[] : nodata;
  boundary({omega_clean});

  stats om_masked = statsf(omega_clean);

  view (fov = 3, camera = "front", tx = -0.05, ty = 0., bg = {1,1,1},
        width = 3000, height = 3000);
  draw_vof ("ibm", "ibmf", filled = -1, lw = 3);
  squares (omega_clean, linear = true, min = om_masked.min, max = om_masked.max);
  save (fname);
}

event adapt (i++) {
  scalar ibmsf[];
  foreach()
    ibmsf[] = vertex_average(point, ibm);
  adapt_wavelet ({ibmsf,u}, (double[]){1.e-15,3e-3,3e-3},
                 maxlevel, minlevel = maxlevel-6);
}

bool is_triple_point (Point point, coord nf, coord ns) {
  return false;
}
bool is_triple_point (Point point, coord nf, coord ns) {
  return false;
}

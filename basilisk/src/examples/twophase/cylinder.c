#include "ibm-gcm.h"
#include "my-centered.h"
#include "my-vof.h"
#include "my-tension.h"
#include "my-two-phase.h"
#include "ibm-gcm-events.h"
#include "view.h"

#define rho_water 999.06
#define mu_water 1.1373e-3
#define rho_air 1.225
#define mu_air 1.789e-5
#define chord 1.0
#define SIGMA 0.0735
#define RATE (sqrt(rho_water*cube(chord)/SIGMA))

double Re = 10000;

double framerate = 1.0/25.0/RATE;
double u0, uref, We;

double film = 0.1;
double t_end = 1.0;
int maxlevel = 11, minlevel = 7;

scalar omega[],m[];

u.n[left] = dirichlet (uref);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (0)
 
event init (t = 0) {
  astats ss;
  int ic = 0;
  do {
    ic++;
    solid (ibm, ibmf, sq(x - chord/2) + sq(y) - sq(chord/2));
    fraction (f, - (sq(x - chord/2) + sq(y) - sq(chord*(1+film)/2)));
    ss = adapt_wavelet ({ibm, f}, (double[]) {1e-30, 1e-30}, maxlevel, minlevel);
  } while ((ss.nf || ss.nc) && ic < 100);
  solid (ibm, ibmf, sq(x - chord/2) + sq(y) - sq(chord/2));
  fraction (f, - (sq(x - chord/2) + sq(y) - sq(chord*(1+film)/2)));
  foreach()
    u.x[] = ibm[];
  boundary((scalar *){u});
}

event output (i++ ; t <= t_end) {
  vorticity (u, omega);
  foreach()
    omega[] = ibm[] < 1. ? nodata : omega[];
  coord Fp, Fmu;
  ibm_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x+Fmu.x)/(0.5*sq((uref))*(chord));
  double CL = (Fp.y+Fmu.y)/(0.5*sq((uref))*(chord));
  stats om  = statsf(omega);
  stats pr  = statsf(p);
  stats ux = statsf(u.x);
  fprintf (stderr, "%d %g %g %g %g %g %g %g\n", i, t, om.max, CL, CD, ux.max, pr.max, pr.min);
  fflush(stderr);
}

event movies (t += framerate, t <= t_end) {
  view(fov = 2, tx = -0.025, width = 1920, height = 1080);
  draw_vof ("f", "fm", filled = 0, lc = {0.0, 0.0, 0.0});
  draw_vof ("ibm", "ibmf", filled = -1, fc = {0.0, 0.0, 0.0});
  stats pr  = statsf(p);
  squares("f", min = 0, max = 1, linear = true);
  save("cylinder.mp4");
}

event adapt (i++) {
  scalar ibmsf[];
  foreach()
    ibmsf[] = vertex_average(point, ibm);
  adapt_wavelet ({ibmsf,f,u}, (double[]){1.e-15,1e-15,3e-3,3e-3}, maxlevel, minlevel);
}

int main (int argc, char *argv[]) {
  if (argc > 1) Re = atoi(argv[1]);
  if (argc > 2) film = atof(argv[2]);
  if (argc > 3) maxlevel = atoi(argv[3]);
  if (argc > 4) t_end    = atof(argv[4]);
  L0 = 16;
  N = 1 << (minlevel);
  origin (-L0/8, -L0/2);
  TOLERANCE = 1e-6;
  f.sigma = SIGMA/SIGMA;
  rho1 = rho_water/rho_water;
  rho2 = rho_air/rho_water;
  mu1 = mu_water/sqrt(rho_water*chord*SIGMA);
  mu2 = mu_air/sqrt(rho_water*chord*SIGMA);
  u0 = Re*mu_air/(rho_air*chord);
  uref = u0/sqrt(SIGMA/(rho_water*chord));
  Re = rho2*uref*chord/mu2;
  We = rho2*sq(uref)*chord/f.sigma;
  if (!pid()) {
    fprintf(stdout, "Re=%.5g\n", Re);
    fprintf(stdout, "We=%.5g\n", We);
    fprintf(stdout, "uref=%.5g\n", uref);
  }
  char log_filename[50];
  snprintf (log_filename, sizeof(log_filename), "cylinder.log");
  freopen (log_filename, "w", stderr);
  run();
}
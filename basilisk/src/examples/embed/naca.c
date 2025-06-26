#include "embed.h"
#include "navier-stokes/centered.h"
#include "trailingcap.h"
#include "view.h"

double mm=0., pp=0., tt=0.12;
double Reynolds = 10000;
double aoa = 2. * M_PI / 180.0;
const double t_end = 30;
const double chord = 1.0;
const double uref = 1.0;
const char* nacaset="0012";
int maxlevel = 13;
int minlevel = 7;
const coord cr = {0.25, 0.0};

face vector muv[];

void nacaset_f(const char* nacaset, double* mm, double* pp, double* tt)
{
    *mm = (nacaset[0] - '0') * 0.01;
    *pp = (nacaset[1] - '0') * 0.1;
    *tt = ((nacaset[2] - '0') * 10 + (nacaset[3] - '0')) * 0.01;
}

double naca(double x, double y, double mm, double pp, double tt) {
  int i = (int)((tt - 0.01) / 0.01);
  double cap = cap_vals[i];
  double te  = r_vals[i];
  double xx = ((x-cr.x)*cos(aoa)-(y-cr.y)*sin(aoa)+cr.x)/chord;
  double yy = ((x-cr.x)*sin(aoa)+(y-cr.y)*cos(aoa)+cr.y)/chord;     
  if (xx >= 0.0 && xx <= 1.0) {     
    double yt = 5.0*tt*(0.2969*sqrt(xx)-0.1260*xx-0.3516*sq(xx)+0.2843*cube(xx)-0.1015*pow(xx, 4.0));
    double yc;
    if (xx >= cap) {
        yt = sqrt(fmax(0.0, sq(te) - sq(xx - (1 - te))));
    }
    if (xx < pp)
      yc = mm/sq(pp)*(2.0*pp*xx - sq(xx));
    else
      yc = mm/sq(1.0 - pp)*(1.0 - 2.0*pp + 2.0*pp*xx - sq(xx));
    return sq(yy - yc) - sq(yt);
    }
  else {
    return 1.0;
  }
}

u.n[left]  = dirichlet(1);
p[left]    = neumann(0);
pf[left]   = neumann(0);

u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/Reynolds;
}

event init (t = 0) {
  scalar csf[];
  astats ss;
  int ic = 0;
  do {
    ic++;
    solid (cs, fs, naca(x, y, mm, pp, tt));
    foreach()
      csf[] = cs[];
    ss = adapt_wavelet ({csf}, (double[]) {1.e-30}, maxlevel, minlevel);
  } while ((ss.nf || ss.nc) && ic < 100);
  solid (cs, fs, naca(x, y, mm, pp, tt));
  foreach()
    u.x[] =  1.;
}

event logfile (i++; t <= t_end) {
  scalar omega[];
  vorticity (u,omega);
  face vector muv[];
  foreach()
    omega[] = cs[] < 1. ? nodata : omega[];
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq((uref))*(chord));
  double CL = (Fp.y + Fmu.y)/(0.5*sq((uref))*(chord));
  stats om  = statsf(omega);
  fprintf (stderr, "%d %g %g %g %g\n", i, t, om.max, CL, CD);
  fflush(stderr);
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-15,3e-3,3e-3}, maxlevel, minlevel);
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
  origin (-L0/8, -L0/2.);
  N = 1 << minlevel;
  TOLERANCE = 1.e-6;
  mu = muv;
  run(); 
}

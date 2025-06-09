#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

double mm=0., pp=0., tt=0.12; // camber,location,thickness
double Reynolds = 10000;
double aoa = 2. * M_PI / 180.0;
const double t_end = 30;

const char* nacaset="0012";
int maxlevel = 13;
face vector muv[];

#define chord   (1.)
#define uref    (1.)
#define tref    ((chord)/(uref))
#define naca00xx(x,y,a) (sq(y) - sq(5.*(a)*(0.2969*sqrt((x)) - 0.1260*((x)) - 0.3516*sq((x)) + 0.2843*cube((x)) - 0.1015*pow((x), 4.))))

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
  origin (-L0/8, -L0/2.);
  N = 1 << (maxlevel-4);
  TOLERANCE = 1.e-6;
  mu = muv;
  
  run(); 
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/Reynolds;
}

u.n[left]  = dirichlet(1);
p[left]    = neumann(0);
pf[left]   = neumann(0);

u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

double naca(double x, double y, double mm, double pp, double tt)
{
  if (x >= 0. && x <= chord) {
    double xr = x * cos(aoa) - y * sin(aoa);
    double yr = x * sin(aoa) + y * cos(aoa);
    double xc = xr / chord;
    double yc = yr / chord;

    if (xc < 0.) xc = 0.;
    if (xc > chord) xc = chord;

    double thetac = 0.;
    if (xc < pp) {
      yc -= mm / sq(pp) * (2. * pp * xc - sq(xc));
      thetac = atan(2. * mm / sq(pp) * (pp - xc));
    } else {
      yc -= mm / sq(1. - pp) * (1. - 2. * pp + 2. * pp * xc - sq(xc));
      thetac = atan(2. * mm / sq(1. - pp) * (pp - xc));
    }

    double geometry = naca00xx(xc, yc, tt * cos(thetac));

    if (xc > 0.98) {
      double d = xc - 0.98;
      geometry += 10. * d * d;
    }

    return geometry;
  }
  else {
    return 1.;
  }
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
    ss = adapt_wavelet ({csf}, (double[]) {1.e-30},
			maxlevel, minlevel = (maxlevel-4));
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
    omega[] = cs[] < 1. ? nodata : fabs(omega[]);
	
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq((uref))*(chord));
  double CL = (Fp.y + Fmu.y)/(0.5*sq((uref))*(chord));
  
  stats om  = statsf(omega);
  fprintf (stderr, "%d %g %g %g %g\n", i, t, om.max, CL, CD);
  fflush(stderr);
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-15,3e-3,3e-3}, maxlevel, 7);
}

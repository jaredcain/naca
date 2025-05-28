// qcc -autolink -Wall -O2 naca.c -o naca -lm -lfb_tiny
// ./naca <NACA> <AoA> <<Re>><<LEVEL>>

#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

double mm=0., pp=0., tt=0.12; // camber,location,thickness
double Reynolds = 10000;
double aoa = 0. * M_PI / 180.0;

const char* nacaset;
int maxlevel = 11;
face vector muv[];

#define chord   (1.)
#define uref    (1.)
#define tref    ((chord)/(uref))
#define naca00xx(x,y,a) (sq(y) - sq(5.*(a)*(0.2969*sqrt((x)) - 0.1260*((x)) - 0.3516*sq((x)) + 0.2843*cube((x)) - 0.1036*pow((x), 4.))))

void nacaset_f(const char* nacaset, double* mm, double* pp, double* tt)
{
    *mm = (nacaset[0] - '0') * 0.01;
    *pp = (nacaset[1] - '0') * 0.1;
    *tt = ((nacaset[2] - '0') * 10 + (nacaset[3] - '0')) * 0.01;
}

int main(int argc, char *argv[])
{
  nacaset = argv[1];
  aoa = atof(argv[2]) * M_PI / 180.0;
  if (argc > 3) {
    Reynolds = atof(argv[3]);
    if (argc > 4) {
      maxlevel = atoi(argv[4]);
    }
  }
  
  nacaset_f(nacaset, &mm, &pp, &tt);
  
  char log_filename[50];
  snprintf(log_filename, sizeof(log_filename), "%s_%.0f.log", nacaset, aoa * 180.0 / M_PI);
  freopen(log_filename, "w", stderr);
  
  L0 = 16.;
  origin (-L0/2, -L0/2.);
  N = 1 << maxlevel;
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
  if (x >= 0. && x <= (chord)) {
    double xr = x * cos(aoa) - y * sin(aoa);
    double yr = x * sin(aoa) + y * cos(aoa);
    double xc = xr / chord;
    double yc = yr / chord;
    if (xc < 0.) {
        xc = 0.;  
    }
    if (xc > chord) {
        xc = chord;  
    }
    double thetac = 0.;
    if (xc < pp) {
      yc -= mm / sq(pp) * (2. * pp * xc - sq(xc));
      thetac = atan(2. * mm / sq(pp) * (pp - xc));
    } else {
      yc -= mm / sq(1. - pp) * (1. - 2. * pp + 2. * pp * xc - sq(xc));
      thetac = atan(2. * mm / sq(1. - pp) * (pp - xc));
    }
    return naca00xx(xc, yc, tt * cos(thetac));
  }
  else {
    return 1.;
  }
}

event init (t = 0)
{
  solid (cs, fs, naca(x,y,mm,pp,tt));
  foreach()
    u.x[] =  1.;
}

event logfile (i += 10){
  scalar omega[];
  vorticity (u,omega);
  face vector muv[];
  foreach()
    {
      if(cs[]< 1.0)
	{
	  omega[]=nodata;
	}
      else
	{
	  omega[]=fabs(omega[]);
	}
    }
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq((uref))*(chord));
  double CL = (Fp.y + Fmu.y)/(0.5*sq((uref))*(chord));
  
  stats om  = statsf(omega);
  fprintf (stderr, "%d %g %g %g %g\n", i, t, om.max, CL, CD);
}

event movies (t += 0.05; t <= 30.)
{
  scalar omega[],m[];
  vorticity(u, omega);
  view(quat = {0.000, 0.000, 0.000, 1.000}, fov = 30, near = 0.01, far = 1000, tx = -0.0, ty = 0.05, tz = -2.25, width = 512, height = 512);
  box();
  squares(color = "omega"); //for some reason freaks out and is a color bomb with out this one
  squares(color = "omega", spread = -1, cbar = true, border = true, pos = {-0.725, 0.4}, mid = true, format = " %.0f", levels = 100, size = 17 , lw=1, fsize = 75);
  //cells();
  draw_vof(c="cs",lw=1, lc = {0,0,0});
  foreach()
    m[] = cs[] - 0.5;
  char zoom[50];
  snprintf(zoom, sizeof(zoom), "%s_%.0f_zoom.mp4", nacaset, aoa * 180.0 / M_PI);
  output_ppm (omega, file = zoom, box = {{-1,-1},{4,1}}, min = -10, max = 10, linear = true, mask = m);
  char video[50];
  snprintf(video, sizeof(video), "%s_%.0f.mp4", nacaset, aoa * 180.0 / M_PI);
  save(video);
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-15,1e-3,1e-3}, maxlevel, 7);
}

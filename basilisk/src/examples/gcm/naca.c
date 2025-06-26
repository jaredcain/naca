#include "ibm-gcm.h"
#include "my-centered.h"
#include "ibm-gcm-events.h"
#include "trailingcap.h"
#include "view.h"

double mm = 0.0, pp = 0.0, tt = 0.12;
double aoa = 2.0*M_PI/180.0;
double Re = 10000;
double t_end = 30;
int maxlevel = 13, minlevel = 7;
const char* nacaset="0012";
const double framerate = 1.0/25.0;
const double chord = 1.0;
const double uref = 1.0;
const coord cr = {0.25, 0.0};

face vector muv[];
scalar omega[],m[],f[];

void nacaset_parse (const char* nacaset, double* mm, double* pp, double* tt) {
  *mm = (nacaset[0]-'0')*0.01;
  *pp = (nacaset[1]-'0')*0.1;
  *tt = ((nacaset[2]-'0')*10+(nacaset[3]-'0'))*0.01;
}

void naca (scalar c, face vector f, double aoa, vertex scalar phii = {0}) {
  int i = (int)((tt - 0.01) / 0.01);
  double cap = cap_vals[i];
  double te  = r_vals[i];
  vertex scalar phi = automatic (phii);
  foreach_vertex() {
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
      phi[] = sq(yy - yc) - sq(yt);
    }
    else {
      phi[] = 1.0;
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

event output (i++ ; t <= t_end) {
  foreach()
    omega[] = ibm[] < 1. ? nodata : fabs(omega[]);
  coord Fp, Fmu;
  ibm_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x+Fmu.x)/(0.5*sq((uref))*(chord));
  double CL = (Fp.y+Fmu.y)/(0.5*sq((uref))*(chord));
  vorticity (u, omega);
  stats om  = statsf(omega);
  stats pr  = statsf(p);
  stats ux = statsf(u.x);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", i, t, om.max, CL, CD, ux.max, pr.max, pr.min);
  fflush(stderr);
}

event movies (t += framerate, t <= t_end) {
  view(fov = 1, tx = -0.025, width = 1920, height = 1080);
  clear();
  draw_vof ("ibm", "ibmf", filled = -1);
  stats pr  = statsf(p);
  squares(color = "p", max = pr.max, min = pr.min); 
  cells();
  char zoom[50];
  snprintf(zoom, sizeof(zoom), "%s_%.0f_%.0f_%d_zoom.mp4", nacaset, aoa * 180.0 / M_PI, Re, maxlevel);
  save(zoom);
  clear();
  view(fov = 21, tx = -0.3375, ty = 0.035, width = 1080, height = 1080);
  box();
  draw_vof ("ibm", "ibmf", filled = -1);
  squares(color = "p", 
	max = pr.max, min = pr.min, 
	cbar = true, 
	border = true, 
	pos = {0.65, 0.6}, 
	label = "p", 
	mid = true, 
	format = " %0.3f", 
	levels = 100, 
	size = 30, 
	lw=1, 
	fsize = 100)
	; 
  char domain[50];
  snprintf(domain, sizeof(domain), "%s_%.0f_%.0f_%d.mp4", nacaset, aoa * 180.0 / M_PI, Re, maxlevel);
  save(domain);
}

event adapt (i++) {
  scalar ibmsf[];
  foreach()
    ibmsf[] = vertex_average(point, ibm);
  adapt_wavelet ({ibmsf,u}, (double[]){1.e-15,3e-3,3e-3}, maxlevel, minlevel);
}

int main (int argc, char *argv[]) {
  if (argc > 1) {
    nacaset = argv[1];
    if (argc > 2) aoa = atof(argv[2])*M_PI/180.0;
    if (argc > 3) Re = atof(argv[3]);
    if (argc > 4) maxlevel = atoi(argv[4]);
	if (argc > 5) t_end = atof(argv[5]);
  }
  nacaset_parse (nacaset, &mm, &pp, &tt);
  char log_filename[50];
  snprintf (log_filename, sizeof(log_filename), "%s_%.0f_%.0f_%d.log", nacaset, aoa * 180.0 / M_PI, Re, maxlevel);
  freopen (log_filename, "w", stderr);
  L0 = 16;
  N = 1 << (minlevel);
  origin (-L0/8, -L0/2);
  TOLERANCE = 1e-6;
  mu = muv;
  run();
}

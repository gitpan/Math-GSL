#include <EXTERN.h>
#include <perl.h>
#include <XSUB.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

MODULE = Math::Gsl::Polynomial  PACKAGE = Math::Gsl::Polynomial	PREFIX = gsl_

PROTOTYPES: DISABLE

void
gsl_poly_complex_solve(...)
PPCODE:
{
 double *a;
 double *z;
 int i;
 int n = items;
 Newz(42,a,items,double);
 for (i=0; i < items; i++)
  {
   a[i] = SvNV(ST(i));
#if 0   
   PerlIO_printf(PerlIO_stderr(),"a[%d]=%g\n",i,a[i]);
#endif   
  }
 while (n && a[n-1] == 0.0)
  n--;
 if (n > 2)
  {  
   gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(n); 
   Newz(42,z,2*(n-1),double);
   i = gsl_poly_complex_solve(a,n,w,z);
   gsl_poly_complex_workspace_free (w);
   Safefree(a);
   if (i == GSL_SUCCESS) 
    {
     for (i=0; i < 2*(n-1); i++)
      {
       XPUSHs(sv_2mortal(newSVnv(z[i])));
#if 0   
       PerlIO_printf(PerlIO_stderr(),"z[%d]=%g\n",i,z[i]);
#endif     
      }
     Safefree(z); 
     XSRETURN(2*(n-1)); 
    }
   else
    {
     Safefree(z); 
    }   
  } 
 XSRETURN(0); 
}




#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"




static int
not_here(char *s)
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

MODULE = Math::Gsl            PACKAGE = Math::Gsl




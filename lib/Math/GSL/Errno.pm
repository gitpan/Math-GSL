# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Errno;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Errnoc;
bootstrap Math::GSL::Errno;
package Math::GSL::Errno;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Errno;

sub TIEHASH {
    my ($classname,$obj) = @_;
    return bless $obj, $classname;
}

sub CLEAR { }

sub FIRSTKEY { }

sub NEXTKEY { }

sub FETCH {
    my ($self,$field) = @_;
    my $member_func = "swig_${field}_get";
    $self->$member_func();
}

sub STORE {
    my ($self,$field,$newval) = @_;
    my $member_func = "swig_${field}_set";
    $self->$member_func($newval);
}

sub this {
    my $ptr = shift;
    return tied(%$ptr);
}


# ------- FUNCTION WRAPPERS --------

package Math::GSL::Errno;

*gsl_error = *Math::GSL::Errnoc::gsl_error;
*gsl_stream_printf = *Math::GSL::Errnoc::gsl_stream_printf;
*gsl_strerror = *Math::GSL::Errnoc::gsl_strerror;
*gsl_set_error_handler = *Math::GSL::Errnoc::gsl_set_error_handler;
*gsl_set_error_handler_off = *Math::GSL::Errnoc::gsl_set_error_handler_off;
*gsl_set_stream_handler = *Math::GSL::Errnoc::gsl_set_stream_handler;
*gsl_set_stream = *Math::GSL::Errnoc::gsl_set_stream;

# ------- VARIABLE STUBS --------

package Math::GSL::Errno;

*GSL_MAJOR_VERSION = *Math::GSL::Errnoc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Errnoc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Errnoc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Errnoc::GSL_NEGZERO;
*GSL_SUCCESS = *Math::GSL::Errnoc::GSL_SUCCESS;
*GSL_FAILURE = *Math::GSL::Errnoc::GSL_FAILURE;
*GSL_CONTINUE = *Math::GSL::Errnoc::GSL_CONTINUE;
*GSL_EDOM = *Math::GSL::Errnoc::GSL_EDOM;
*GSL_ERANGE = *Math::GSL::Errnoc::GSL_ERANGE;
*GSL_EFAULT = *Math::GSL::Errnoc::GSL_EFAULT;
*GSL_EINVAL = *Math::GSL::Errnoc::GSL_EINVAL;
*GSL_EFAILED = *Math::GSL::Errnoc::GSL_EFAILED;
*GSL_EFACTOR = *Math::GSL::Errnoc::GSL_EFACTOR;
*GSL_ESANITY = *Math::GSL::Errnoc::GSL_ESANITY;
*GSL_ENOMEM = *Math::GSL::Errnoc::GSL_ENOMEM;
*GSL_EBADFUNC = *Math::GSL::Errnoc::GSL_EBADFUNC;
*GSL_ERUNAWAY = *Math::GSL::Errnoc::GSL_ERUNAWAY;
*GSL_EMAXITER = *Math::GSL::Errnoc::GSL_EMAXITER;
*GSL_EZERODIV = *Math::GSL::Errnoc::GSL_EZERODIV;
*GSL_EBADTOL = *Math::GSL::Errnoc::GSL_EBADTOL;
*GSL_ETOL = *Math::GSL::Errnoc::GSL_ETOL;
*GSL_EUNDRFLW = *Math::GSL::Errnoc::GSL_EUNDRFLW;
*GSL_EOVRFLW = *Math::GSL::Errnoc::GSL_EOVRFLW;
*GSL_ELOSS = *Math::GSL::Errnoc::GSL_ELOSS;
*GSL_EROUND = *Math::GSL::Errnoc::GSL_EROUND;
*GSL_EBADLEN = *Math::GSL::Errnoc::GSL_EBADLEN;
*GSL_ENOTSQR = *Math::GSL::Errnoc::GSL_ENOTSQR;
*GSL_ESING = *Math::GSL::Errnoc::GSL_ESING;
*GSL_EDIVERGE = *Math::GSL::Errnoc::GSL_EDIVERGE;
*GSL_EUNSUP = *Math::GSL::Errnoc::GSL_EUNSUP;
*GSL_EUNIMPL = *Math::GSL::Errnoc::GSL_EUNIMPL;
*GSL_ECACHE = *Math::GSL::Errnoc::GSL_ECACHE;
*GSL_ETABLE = *Math::GSL::Errnoc::GSL_ETABLE;
*GSL_ENOPROG = *Math::GSL::Errnoc::GSL_ENOPROG;
*GSL_ENOPROGJ = *Math::GSL::Errnoc::GSL_ENOPROGJ;
*GSL_ETOLF = *Math::GSL::Errnoc::GSL_ETOLF;
*GSL_ETOLX = *Math::GSL::Errnoc::GSL_ETOLX;
*GSL_ETOLG = *Math::GSL::Errnoc::GSL_ETOLG;
*GSL_EOF = *Math::GSL::Errnoc::GSL_EOF;
 
@EXPORT_OK = qw/
                gsl_error gsl_stream_printf gsl_strerror gsl_set_error_handler 
                gsl_set_error_handler_off gsl_set_stream_handler gsl_set_stream 
                $GSL_SUCCESS $GSL_FAILURE $GSL_CONTINUE $GSL_EDOM 
                $GSL_ERANGE $GSL_EFAULT $GSL_EINVAL $GSL_EFAILED 
                $GSL_EFACTOR $GSL_ESANITY $GSL_ENOMEM $GSL_EBADFUNC 
                $GSL_ERUNAWAY $GSL_EMAXITER $GSL_EZERODIV $GSL_EBADTOL 
                $GSL_ETOL $GSL_EUNDRFLW $GSL_EOVRFLW $GSL_ELOSS 
                $GSL_EROUND $GSL_EBADLEN $GSL_ENOTSQR $GSL_ESING 
                $GSL_EDIVERGE $GSL_EUNSUP $GSL_EUNIMPL $GSL_ECACHE 
                $GSL_ETABLE $GSL_ENOPROG $GSL_ENOPROGJ $GSL_ETOLF 
                $GSL_ETOLX $GSL_ETOLG $GSL_EOF 
            /;
%EXPORT_TAGS = ( all => \@EXPORT_OK );

=head1 NAME

Math::GSL::Errno - Error Handling

=cut

=head1 SYNOPSIS

    use Math::GSL::Errno qw/:all/;
    gsl_set_error_handler_off();

    print gsl_strerror($GSL_EDOM) . "\n";

=head1 STATUS CODES

=over 4

=item * $GSL_SUCCESS 

Success

=item * $GSL_FAILURE 

General Failure

=item * $GSL_CONTINUE 

Iteration has not converged

=item * $GSL_EDOM

Domain error; used by mathematical functions when an argument value does not fall into the domain over which the function is defined (like EDOM in the C library)

=item * $GSL_ERANGE

Range error; used by mathematical functions when the result value is not representable because of overflow or underflow (like ERANGE in the C library) 

=item * $GSL_EFAULT

Invalid Pointer

=item * $GSL_EINVAL

Invalid argument. This is used to indicate various kinds of problems with passing the wrong argument to a library function (like EINVAL in the C library).Invalid argument. This is used to indicate various kinds of problems with passing the wrong argument to a library function (like EINVAL in the C library). 

=item * $GSL_EFAILED

Generic Failure

=item * $GSL_EFACTOR

Factorization Failed

=item * $GSL_ESANITY

Sanity Check Failed

=item * $GSL_ENOMEM

No memory available. The system cannot allocate more virtual memory because its capacity is full (like ENOMEM in the C library). This error is reported when a GSL routine encounters problems when trying to allocate memory with malloc.

=item * $GSL_EBADFUNC

Problem with user-supplied function

=item * $GSL_ERUNAWAY

Iterative process is our of control

=item * $GSL_EMAXITER

Exceeded max number of iterations

=item * $GSL_EZERODIV

Division by zero

=item * $GSL_EBADTOL

Invalid user-specified tolerance

=item * $GSL_ETOL

Failed to reach the specified tolerance

=item * $GSL_EUNDRFLW

Underflow

=item * $GSL_EOVRFLW

Overflow

=item * $GSL_ELOSS

Loss of accuracy

=item * $GSL_EROUND

Failed due to roundoff error

=item * $GSL_EBADLEN

Matrix/vector lengths not compatible

=item * $GSL_ENOTSQR

Not a square matrix

=item * $GSL_ESING

Singularity Detected

=item * $GSL_EDIVERGE

Integral/Series is divergent

=item * $GSL_EUNSUP

Not supported by hardware

=item * $GSL_EUNIMPL

Not implemented

=item * $GSL_ECACHE

Cache limit exceeded

=item * $GSL_ETABLE

Table limit exceeded

=item * $GSL_ENOPROG

Iteration not converging

=item * $GSL_ENOPROGJ

Jacobian not improving solution

=item * $GSL_ETOLF

Cannot reach tolerance in F

=item * $GSL_ETOLX

Cannot reach tolerance in X

=item * $GSL_ETOLG

Cannot reach tolerance in Gradient

=item * $GSL_EOF

End of file

=back

=head1 FUNCTIONS

=over

=item * gsl_error

=item * gsl_stream_printf

=item * gsl_strerror($gsl_errno) - This function returns a pointer to a string describing the error code gsl_errno. For example, print ("error: gsl_strerror ($status)\n"); would print an error message like error: output range error for a status value of GSL_ERANGE.

=item * gsl_set_error_handler 

=item * gsl_set_error_handler_off() - This function turns off the error handler by defining an error handler which does nothing. This will cause the program to continue after any error, so the return values from any library routines must be checked. This is the recommended behavior for production programs. The previous handler is returned (so that you can restore it later).

=item * gsl_set_stream_handler

=item * gsl_set_stream 

=back

=cut

1;

1;

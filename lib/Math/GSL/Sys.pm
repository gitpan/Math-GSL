# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Sys;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Sysc;
bootstrap Math::GSL::Sys;
package Math::GSL::Sys;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Sys;

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

package Math::GSL::Sys;

*gsl_log1p = *Math::GSL::Sysc::gsl_log1p;
*gsl_expm1 = *Math::GSL::Sysc::gsl_expm1;
*gsl_hypot = *Math::GSL::Sysc::gsl_hypot;
*gsl_hypot3 = *Math::GSL::Sysc::gsl_hypot3;
*gsl_acosh = *Math::GSL::Sysc::gsl_acosh;
*gsl_asinh = *Math::GSL::Sysc::gsl_asinh;
*gsl_atanh = *Math::GSL::Sysc::gsl_atanh;
*gsl_isnan = *Math::GSL::Sysc::gsl_isnan;
*gsl_isinf = *Math::GSL::Sysc::gsl_isinf;
*gsl_finite = *Math::GSL::Sysc::gsl_finite;
*gsl_nan = *Math::GSL::Sysc::gsl_nan;
*gsl_posinf = *Math::GSL::Sysc::gsl_posinf;
*gsl_neginf = *Math::GSL::Sysc::gsl_neginf;
*gsl_fdiv = *Math::GSL::Sysc::gsl_fdiv;
*gsl_coerce_double = *Math::GSL::Sysc::gsl_coerce_double;
*gsl_coerce_float = *Math::GSL::Sysc::gsl_coerce_float;
*gsl_coerce_long_double = *Math::GSL::Sysc::gsl_coerce_long_double;
*gsl_ldexp = *Math::GSL::Sysc::gsl_ldexp;
*gsl_frexp = *Math::GSL::Sysc::gsl_frexp;
*gsl_fcmp = *Math::GSL::Sysc::gsl_fcmp;

# ------- VARIABLE STUBS --------

package Math::GSL::Sys;

*GSL_MAJOR_VERSION = *Math::GSL::Sysc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Sysc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Sysc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Sysc::GSL_NEGZERO;

our @EXPORT = qw();
our @EXPORT_OK = qw/
               gsl_log1p
               gsl_expm1
               gsl_hypot
               gsl_hypot3
               gsl_acosh
               gsl_asinh
               gsl_atanh
               gsl_isnan
               gsl_isinf
               gsl_finite
               gsl_posinf
               gsl_neginf
               gsl_fdiv
               gsl_coerce_double
               gsl_coerce_float
               gsl_coerce_long_double
               gsl_ldexp
               gsl_frexp
               gsl_fcmp
               gsl_nan
               gsl_isnan
               gsl_inf
               $GSL_NAN
               $GSL_POSINF
               $GSL_NEGINF
             /;

our %EXPORT_TAGS = ( all => \@EXPORT_OK );
our $GSL_NAN    = gsl_nan();
our $GSL_POSINF = gsl_posinf();
our $GSL_NEGINF = gsl_neginf();

__END__

=head1 NAME

Math::GSL::Sys - Misc Math Functions

=head1 SYNOPSIS

    use Math::GSL::Sys qw/:all/;

=head1 DESCRIPTION

This module contains various useful math functions that are not usually
provided by standard libraries.

=over

=item * C<gsl_log1p($x)> 

This function computes the value of \log(1+$x) in a way that is accurate for
small $x. It provides an alternative to the BSD math function log1p(x).

=item * C<gsl_expm1($x)> 

This function computes the value of \exp($x)-1 in a way that is accurate for
small $x. It provides an alternative to the BSD math function expm1(x).

=item * C<gsl_hypot($x, $y)> 

This function computes the value of \sqrt{$x^2 + $y^2} in a way that avoids
overflow. It provides an alternative to the BSD math function hypot($x,$y).

=item * C<gsl_hypot3($x, $y, $z)> 

This function computes the value of \sqrt{$x^2 + $y^2 + $z^2} in a way that
avoids overflow.

=item * C<gsl_acosh($x)> 

This function computes the value of \arccosh($x). It provides an alternative to
the standard math function acosh($x).

=item * C<gsl_asinh($x)> 

This function computes the value of \arcsinh($x). It provides an alternative to
the standard math function asinh($x).

=item * C<gsl_atanh($x)> 

This function computes the value of \arctanh($x). It provides an alternative to
the standard math function atanh($x).

=item * C<gsl_isnan($x)> 

This function returns 1 if $x is not-a-number.

=item * C<gsl_isinf($x)>  

This function returns +1 if $x is positive infinity, -1 if $x is negative
infinity and 0 otherwise.

=item * C<gsl_finite($x)> 

This function returns 1 if $x is a real number, and 0 if it is infinite or not-a-number.

=item * C<gsl_posinf >

=item * C<gsl_neginf >

=item * C<gsl_fdiv >

=item * C<gsl_coerce_double >

=item * C<gsl_coerce_float >

=item * C<gsl_coerce_long_double >

=item * C<gsl_ldexp($x, $e)> 

This function computes the value of $x * 2**$e. It provides an alternative to
the standard math function ldexp($x,$e).

=item * C<gsl_frexp($x)> 

This function splits the number $x into its normalized fraction f and exponent
e, such that $x = f * 2^e and 0.5 <= f < 1. The function returns f and then the
exponent in e. If $x is zero, both f and e are set to zero. This function
provides an alternative to the standard math function frexp(x, e).

=item * C<gsl_fcmp($x, $y, $epsilon)> 

This function determines whether $x and $y are approximately equal to a
relative accuracy $epsilon. The relative accuracy is measured using an interval
of size 2 \delta, where \delta = 2^k \epsilon and k is the maximum base-2
exponent of $x and $y as computed by the function frexp. If $x and $y lie
within this interval, they are considered approximately equal and the function
returns 0. Otherwise if $x < $y, the function returns -1, or if $x > $y, the
function returns +1. Note that $x and $y are compared to relative accuracy, so
this function is not suitable for testing whether a value is approximately
zero. The implementation is based on the package fcmp by T.C. Belding.

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

=head1 AUTHORS

Jonathan "Duke" Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2011 Jonathan "Duke" Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

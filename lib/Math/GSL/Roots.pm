# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Roots;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Rootsc;
bootstrap Math::GSL::Roots;
package Math::GSL::Roots;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Roots;

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

package Math::GSL::Roots;

*gsl_root_fsolver_alloc = *Math::GSL::Rootsc::gsl_root_fsolver_alloc;
*gsl_root_fsolver_free = *Math::GSL::Rootsc::gsl_root_fsolver_free;
*gsl_root_fsolver_set = *Math::GSL::Rootsc::gsl_root_fsolver_set;
*gsl_root_fsolver_iterate = *Math::GSL::Rootsc::gsl_root_fsolver_iterate;
*gsl_root_fsolver_name = *Math::GSL::Rootsc::gsl_root_fsolver_name;
*gsl_root_fsolver_root = *Math::GSL::Rootsc::gsl_root_fsolver_root;
*gsl_root_fsolver_x_lower = *Math::GSL::Rootsc::gsl_root_fsolver_x_lower;
*gsl_root_fsolver_x_upper = *Math::GSL::Rootsc::gsl_root_fsolver_x_upper;
*gsl_root_fdfsolver_alloc = *Math::GSL::Rootsc::gsl_root_fdfsolver_alloc;
*gsl_root_fdfsolver_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_set;
*gsl_root_fdfsolver_iterate = *Math::GSL::Rootsc::gsl_root_fdfsolver_iterate;
*gsl_root_fdfsolver_free = *Math::GSL::Rootsc::gsl_root_fdfsolver_free;
*gsl_root_fdfsolver_name = *Math::GSL::Rootsc::gsl_root_fdfsolver_name;
*gsl_root_fdfsolver_root = *Math::GSL::Rootsc::gsl_root_fdfsolver_root;
*gsl_root_test_interval = *Math::GSL::Rootsc::gsl_root_test_interval;
*gsl_root_test_residual = *Math::GSL::Rootsc::gsl_root_test_residual;
*gsl_root_test_delta = *Math::GSL::Rootsc::gsl_root_test_delta;

############# Class : Math::GSL::Roots::gsl_root_fsolver_type ##############

package Math::GSL::Roots::gsl_root_fsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_name_get;
*swig_name_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_name_set;
*swig_size_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_size_get;
*swig_size_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_size_set;
*swig_set_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_set_get;
*swig_set_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_iterate_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fsolver_type($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Math::GSL::Roots::gsl_root_fsolver ##############

package Math::GSL::Roots::gsl_root_fsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_get;
*swig_type_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_set;
*swig_function_get = *Math::GSL::Rootsc::gsl_root_fsolver_function_get;
*swig_function_set = *Math::GSL::Rootsc::gsl_root_fsolver_function_set;
*swig_root_get = *Math::GSL::Rootsc::gsl_root_fsolver_root_get;
*swig_root_set = *Math::GSL::Rootsc::gsl_root_fsolver_root_set;
*swig_x_lower_get = *Math::GSL::Rootsc::gsl_root_fsolver_x_lower_get;
*swig_x_lower_set = *Math::GSL::Rootsc::gsl_root_fsolver_x_lower_set;
*swig_x_upper_get = *Math::GSL::Rootsc::gsl_root_fsolver_x_upper_get;
*swig_x_upper_set = *Math::GSL::Rootsc::gsl_root_fsolver_x_upper_set;
*swig_state_get = *Math::GSL::Rootsc::gsl_root_fsolver_state_get;
*swig_state_set = *Math::GSL::Rootsc::gsl_root_fsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fsolver($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Math::GSL::Roots::gsl_root_fdfsolver_type ##############

package Math::GSL::Roots::gsl_root_fdfsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_name_get;
*swig_name_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_name_set;
*swig_size_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_size_get;
*swig_size_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_size_set;
*swig_set_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_set_get;
*swig_set_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_iterate_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fdfsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fdfsolver_type($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Math::GSL::Roots::gsl_root_fdfsolver ##############

package Math::GSL::Roots::gsl_root_fdfsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_get;
*swig_type_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_set;
*swig_fdf_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_fdf_get;
*swig_fdf_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_fdf_set;
*swig_root_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_root_get;
*swig_root_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_root_set;
*swig_state_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_state_get;
*swig_state_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fdfsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fdfsolver($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


# ------- VARIABLE STUBS --------

package Math::GSL::Roots;

*GSL_MAJOR_VERSION = *Math::GSL::Rootsc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Rootsc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Rootsc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Rootsc::GSL_NEGZERO;

my %__gsl_root_fsolver_bisection_hash;
tie %__gsl_root_fsolver_bisection_hash,"Math::GSL::Roots::gsl_root_fsolver_type", $Math::GSL::Rootsc::gsl_root_fsolver_bisection;
$gsl_root_fsolver_bisection= \%__gsl_root_fsolver_bisection_hash;
bless $gsl_root_fsolver_bisection, Math::GSL::Roots::gsl_root_fsolver_type;

my %__gsl_root_fsolver_brent_hash;
tie %__gsl_root_fsolver_brent_hash,"Math::GSL::Roots::gsl_root_fsolver_type", $Math::GSL::Rootsc::gsl_root_fsolver_brent;
$gsl_root_fsolver_brent= \%__gsl_root_fsolver_brent_hash;
bless $gsl_root_fsolver_brent, Math::GSL::Roots::gsl_root_fsolver_type;

my %__gsl_root_fsolver_falsepos_hash;
tie %__gsl_root_fsolver_falsepos_hash,"Math::GSL::Roots::gsl_root_fsolver_type", $Math::GSL::Rootsc::gsl_root_fsolver_falsepos;
$gsl_root_fsolver_falsepos= \%__gsl_root_fsolver_falsepos_hash;
bless $gsl_root_fsolver_falsepos, Math::GSL::Roots::gsl_root_fsolver_type;

my %__gsl_root_fdfsolver_newton_hash;
tie %__gsl_root_fdfsolver_newton_hash,"Math::GSL::Roots::gsl_root_fdfsolver_type", $Math::GSL::Rootsc::gsl_root_fdfsolver_newton;
$gsl_root_fdfsolver_newton= \%__gsl_root_fdfsolver_newton_hash;
bless $gsl_root_fdfsolver_newton, Math::GSL::Roots::gsl_root_fdfsolver_type;

my %__gsl_root_fdfsolver_secant_hash;
tie %__gsl_root_fdfsolver_secant_hash,"Math::GSL::Roots::gsl_root_fdfsolver_type", $Math::GSL::Rootsc::gsl_root_fdfsolver_secant;
$gsl_root_fdfsolver_secant= \%__gsl_root_fdfsolver_secant_hash;
bless $gsl_root_fdfsolver_secant, Math::GSL::Roots::gsl_root_fdfsolver_type;

my %__gsl_root_fdfsolver_steffenson_hash;
tie %__gsl_root_fdfsolver_steffenson_hash,"Math::GSL::Roots::gsl_root_fdfsolver_type", $Math::GSL::Rootsc::gsl_root_fdfsolver_steffenson;
$gsl_root_fdfsolver_steffenson= \%__gsl_root_fdfsolver_steffenson_hash;
bless $gsl_root_fdfsolver_steffenson, Math::GSL::Roots::gsl_root_fdfsolver_type;

@EXPORT_OK = qw/
               gsl_root_fsolver_alloc 
               gsl_root_fsolver_free 
               gsl_root_fsolver_set 
               gsl_root_fsolver_iterate 
               gsl_root_fsolver_name 
               gsl_root_fsolver_root 
               gsl_root_fsolver_x_lower 
               gsl_root_fsolver_x_upper 
               gsl_root_fdfsolver_alloc 
               gsl_root_fdfsolver_set 
               gsl_root_fdfsolver_iterate 
               gsl_root_fdfsolver_free 
               gsl_root_fdfsolver_name 
               gsl_root_fdfsolver_root 
               gsl_root_test_interval 
               gsl_root_test_residual 
               gsl_root_test_delta 
               $gsl_root_fsolver_bisection    
               $gsl_root_fsolver_brent   
               $gsl_root_fsolver_falsepos     
               $gsl_root_fdfsolver_newton     
               $gsl_root_fdfsolver_secant     
               $gsl_root_fdfsolver_steffenson 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Roots - Routines for finding roots of arbitrary one-dimensional functions.

=head1 SYNOPSIS

use Math::GSL::Roots qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over 

=item * C<gsl_root_fsolver_alloc($T)> - This function returns a pointer to a newly allocated instance of a solver of type $T. $T must be one of the constant included with this module. If there is insufficient memory to create the solver then the function returns a null pointer and the error handler is invoked with an error code of $GSL_ENOMEM.

=item * C<gsl_root_fsolver_free($s)> - This function frees all the memory associated with the solver $s. 

=item * C<gsl_root_fsolver_set($s, $f, $x_lower, $x_upper)> - This function initializes, or reinitializes, an existing solver $s to use the function $f and the initial search interval [$x_lower, $x_upper]. $f has to be of this form : sub { my $x=shift; function_with_$x }. For example, sub { my $x=shift; ($x-3.2)**3 } is a valid value for $f.

=item * C<gsl_root_fsolver_iterate($s)> - This function performs a single iteration of the solver $s. If the iteration encounters an unexpected problem then an error code will be returned (the Math::GSL::Errno has to be included),
 $GSL_EBADFUNC - the iteration encountered a singular point where the function or its derivative evaluated to Inf or NaN.
 $GSL_EZERODIV - the derivative of the function vanished at the iteration point, preventing the algorithm from continuing without a division by zero. 

=item * C<gsl_root_fsolver_name($s)> - This function returns the name of the solver use within the $s solver.

=item * C<gsl_root_fsolver_root($s)> - This function returns the current estimate of the root for the solver $s.

=item * C<gsl_root_fsolver_x_lower($s)> - This function returns the current lower value of the bracketing interval for the solver $s. 

=item * C<gsl_root_fsolver_x_upper($s)> - This function returns the current lower value of the bracketing interval for the solver $s.

=item * C<gsl_root_fdfsolver_alloc($T)> - This function returns a pointer to a newly allocated instance of a derivative-based solver of type $T. If there is insufficient memory to create the solver then the function returns a null pointer and the error handler is invoked with an error code of $GSL_ENOMEM.

=item * C<gsl_root_fdfsolver_set($s, $fdf, $root)> - This function initializes, or reinitializes, an existing solver $s to use the function and derivative $fdf and the initial guess $root. $f has to be of this form : sub { my $x=shift; function_with_$x }. For example, sub { my $x=shift; ($x-3.2)**3 } is a valid value for $fdf.

=item * C<gsl_root_fdfsolver_iterate($s)> - This function performs a single iteration of the solver $s. If the iteration encounters an unexpected problem then an error code will be returned (the Math::GSL::Errno has to be included),
 $GSL_EBADFUNC - the iteration encountered a singular point where the function or its derivative evaluated to Inf or NaN.
 $GSL_EZERODIV - the derivative of the function vanished at the iteration point, preventing the algorithm from continuing without a division by zero. 

=item * C<gsl_root_fdfsolver_free($s)> - This function frees all the memory associated with the solver $s.

=item * C<gsl_root_fdfsolver_name($s)> - This function returns the name of the solver use within the $s solver.

=item * C<gsl_root_fdfsolver_root($s)> - This function returns the current estimate of the root for the solver $s.

=item * C<gsl_root_test_interval($x_lower, $x_upper, $epsabs, $epsrel)> - This function tests for the convergence of the interval [$x_lower, $x_upper] with absolute error epsabs and relative error $epsrel. The test returns $GSL_SUCCESS if the following condition is achieved,
    		|a - b| < epsabs + epsrel min(|a|,|b|)
 when the interval x = [a,b] does not include the origin. If the interval includes the origin then \min(|a|,|b|) is replaced by zero (which is the minimum value of |x| over the interval). This ensures that the relative error is accurately estimated for roots close to the origin.
 This condition on the interval also implies that any estimate of the root r in the interval satisfies the same condition with respect to the true root r^*,
              |r - r^*| < epsabs + epsrel r^*
  assuming that the true root r^* is contained within the interval. 

=item * C<gsl_root_test_residual($f, $epsabs)> - This function tests the residual value $f against the absolute error bound $epsabs. The test returns $GSL_SUCCESS if the following condition is achieved,
              |$f| < $epsabs
    and returns $GSL_CONTINUE otherwise. This criterion is suitable for situations where the precise location of the root, x, is unimportant provided a value can be found where the residual, |f(x)|, is small enough. 

=item * C<gsl_root_test_delta($x1, $x0, $epsabs, $epsrel)> - This function tests for the convergence of the sequence ..., $x0, $x1 with absolute error $epsabs and relative error $epsrel. The test returns $GSL_SUCCESS if the following condition is achieved,
              |x_1 - x_0| < epsabs + epsrel |x_1|
    and returns $GSL_CONTINUE otherwise. 

=back

This module also includes the following constants :

=over

=item * C<$gsl_root_fsolver_bisection> - The bisection algorithm is the simplest method of bracketing the roots of a function. It is the slowest algorithm provided by the library, with linear convergence. On each iteration, the interval is bisected and the value of the function at the midpoint is calculated. The sign of this value is used to determine which half of the interval does not contain a root. That half is discarded to give a new, smaller interval containing the root. This procedure can be continued indefinitely until the interval is sufficiently small. At any time the current estimate of the root is taken as the midpoint of the interval. 

=item * C<$gsl_root_fsolver_brent> - The Brent-Dekker method (referred to here as Brent's method) combines an interpolation strategy with the bisection algorithm. This produces a fast algorithm which is still robust. On each iteration Brent's method approximates the function using an interpolating curve. On the first iteration this is a linear interpolation of the two endpoints. For subsequent iterations the algorithm uses an inverse quadratic fit to the last three points, for higher accuracy. The intercept of the interpolating curve with the x-axis is taken as a guess for the root. If it lies within the bounds of the current interval then the interpolating point is accepted, and used to generate a smaller interval. If the interpolating point is not accepted then the algorithm falls back to an ordinary bisection step. The best estimate of the root is taken from the most recent interpolation or bisection. 

=item * C<$gsl_root_fsolver_falsepos> - The false position algorithm is a method of finding roots based on linear interpolation. Its convergence is linear, but it is usually faster than bisection. On each iteration a line is drawn between the endpoints (a,f(a)) and (b,f(b)) and the point where this line crosses the x-axis taken as a “midpoint”. The value of the function at this point is calculated and its sign is used to determine which side of the interval does not contain a root. That side is discarded to give a new, smaller interval containing the root. This procedure can be continued indefinitely until the interval is sufficiently small. The best estimate of the root is taken from the linear interpolation of the interval on the current iteration. 

=item * C<$gsl_root_fdfsolver_newton> - Newton's Method is the standard root-polishing algorithm. The algorithm begins with an initial guess for the location of the root. On each iteration, a line tangent to the function f is drawn at that position. The point where this line crosses the x-axis becomes the new guess. The iteration is defined by the following sequence, x_{i+1} = x_i - f(x_i)/f'(x_i) Newton's method converges quadratically for single roots, and linearly for multiple roots. 

=item * C<$gsl_root_fdfsolver_secant> - The secant method is a simplified version of Newton's method which does not require the computation of the derivative on every step.
On its first iteration the algorithm begins with Newton's method, using the derivative to compute a first step,
          x_1 = x_0 - f(x_0)/f'(x_0)
 Subsequent iterations avoid the evaluation of the derivative by replacing it with a numerical estimate, the slope of the line through the previous two points,
          x_{i+1} = x_i f(x_i) / f'_{est} where
           f'_{est} = (f(x_i) - f(x_{i-1})/(x_i - x_{i-1})
 When the derivative does not change significantly in the vicinity of the root the secant method gives a useful saving. Asymptotically the secant method is faster than Newton's method whenever the cost of evaluating the derivative is more than 0.44 times the cost of evaluating the function itself. As with all methods of computing a numerical derivative the estimate can suffer from cancellation errors if the separation of the points becomes too small.

 On single roots, the method has a convergence of order (1 + \sqrt 5)/2 (approximately 1.62). It converges linearly for multiple roots. 

=item * C<$gsl_root_fdfsolver_steffenson> - The Steffenson Method provides the fastest convergence of all the routines. It combines the basic Newton algorithm with an Aitken “delta-squared” acceleration. If the Newton iterates are x_i then the acceleration procedure generates a new sequence R_i,
    R_i = x_i - (x_{i+1} - x_i)^2 / (x_{i+2} - 2 x_{i+1} + x_{i})
 which converges faster than the original sequence under reasonable conditions. The new sequence requires three terms before it can produce its first value so the method returns accelerated values on the second and subsequent iterations. On the first iteration it returns the ordinary Newton estimate. The Newton iterate is also returned if the denominator of the acceleration term ever becomes zero.

As with all acceleration procedures this method can become unstable if the function is not well-behaved. 

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

 Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want


=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

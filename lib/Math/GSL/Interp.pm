# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Interp;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Interpc;
bootstrap Math::GSL::Interp;
package Math::GSL::Interp;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Interp;

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

package Math::GSL::Interp;

*gsl_interp_accel_alloc = *Math::GSL::Interpc::gsl_interp_accel_alloc;
*gsl_interp_accel_reset = *Math::GSL::Interpc::gsl_interp_accel_reset;
*gsl_interp_accel_free = *Math::GSL::Interpc::gsl_interp_accel_free;
*gsl_interp_alloc = *Math::GSL::Interpc::gsl_interp_alloc;
*gsl_interp_init = *Math::GSL::Interpc::gsl_interp_init;
*gsl_interp_name = *Math::GSL::Interpc::gsl_interp_name;
*gsl_interp_min_size = *Math::GSL::Interpc::gsl_interp_min_size;
*gsl_interp_eval_e = *Math::GSL::Interpc::gsl_interp_eval_e;
*gsl_interp_eval = *Math::GSL::Interpc::gsl_interp_eval;
*gsl_interp_eval_deriv_e = *Math::GSL::Interpc::gsl_interp_eval_deriv_e;
*gsl_interp_eval_deriv = *Math::GSL::Interpc::gsl_interp_eval_deriv;
*gsl_interp_eval_deriv2_e = *Math::GSL::Interpc::gsl_interp_eval_deriv2_e;
*gsl_interp_eval_deriv2 = *Math::GSL::Interpc::gsl_interp_eval_deriv2;
*gsl_interp_eval_integ_e = *Math::GSL::Interpc::gsl_interp_eval_integ_e;
*gsl_interp_eval_integ = *Math::GSL::Interpc::gsl_interp_eval_integ;
*gsl_interp_free = *Math::GSL::Interpc::gsl_interp_free;
*gsl_interp_bsearch = *Math::GSL::Interpc::gsl_interp_bsearch;
*gsl_interp_accel_find = *Math::GSL::Interpc::gsl_interp_accel_find;

############# Class : Math::GSL::Interp::gsl_interp_accel ##############

package Math::GSL::Interp::gsl_interp_accel;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Interp );
%OWNER = ();
%ITERATORS = ();
*swig_cache_get = *Math::GSL::Interpc::gsl_interp_accel_cache_get;
*swig_cache_set = *Math::GSL::Interpc::gsl_interp_accel_cache_set;
*swig_miss_count_get = *Math::GSL::Interpc::gsl_interp_accel_miss_count_get;
*swig_miss_count_set = *Math::GSL::Interpc::gsl_interp_accel_miss_count_set;
*swig_hit_count_get = *Math::GSL::Interpc::gsl_interp_accel_hit_count_get;
*swig_hit_count_set = *Math::GSL::Interpc::gsl_interp_accel_hit_count_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Interpc::new_gsl_interp_accel(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Interpc::delete_gsl_interp_accel($self);
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


############# Class : Math::GSL::Interp::gsl_interp ##############

package Math::GSL::Interp::gsl_interp;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Interp );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Interpc::gsl_interp_type_get;
*swig_type_set = *Math::GSL::Interpc::gsl_interp_type_set;
*swig_xmin_get = *Math::GSL::Interpc::gsl_interp_xmin_get;
*swig_xmin_set = *Math::GSL::Interpc::gsl_interp_xmin_set;
*swig_xmax_get = *Math::GSL::Interpc::gsl_interp_xmax_get;
*swig_xmax_set = *Math::GSL::Interpc::gsl_interp_xmax_set;
*swig_size_get = *Math::GSL::Interpc::gsl_interp_size_get;
*swig_size_set = *Math::GSL::Interpc::gsl_interp_size_set;
*swig_state_get = *Math::GSL::Interpc::gsl_interp_state_get;
*swig_state_set = *Math::GSL::Interpc::gsl_interp_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Interpc::new_gsl_interp(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Interpc::delete_gsl_interp($self);
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

package Math::GSL::Interp;

*GSL_MAJOR_VERSION = *Math::GSL::Interpc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Interpc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Interpc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Interpc::GSL_NEGZERO;
*gsl_interp_linear = *Math::GSL::Interpc::gsl_interp_linear;
*gsl_interp_polynomial = *Math::GSL::Interpc::gsl_interp_polynomial;
*gsl_interp_cspline = *Math::GSL::Interpc::gsl_interp_cspline;
*gsl_interp_cspline_periodic = *Math::GSL::Interpc::gsl_interp_cspline_periodic;
*gsl_interp_akima = *Math::GSL::Interpc::gsl_interp_akima;
*gsl_interp_akima_periodic = *Math::GSL::Interpc::gsl_interp_akima_periodic;

@EXPORT_OK = qw/
               gsl_interp_accel_alloc 
               gsl_interp_accel_find 
               gsl_interp_accel_reset 
               gsl_interp_accel_free 
               gsl_interp_alloc 
               gsl_interp_init 
               gsl_interp_name 
               gsl_interp_min_size 
               gsl_interp_eval_e 
               gsl_interp_eval 
               gsl_interp_eval_deriv_e 
               gsl_interp_eval_deriv 
               gsl_interp_eval_deriv2_e 
               gsl_interp_eval_deriv2 
               gsl_interp_eval_integ_e 
               gsl_interp_eval_integ 
               gsl_interp_free 
               gsl_interp_bsearch
               $gsl_interp_linear
               $gsl_interp_polynomial
               $gsl_interp_cspline
               $gsl_interp_cspline_periodic
               $gsl_interp_akima
               $gsl_interp_akima_periodic
             /;
%EXPORT_TAGS = ( all => \@EXPORT_OK  );

__END__

=head1 NAME

Math::GSL::Interp - Interpolation

=head1 SYNOPSIS

 use Math::GSL::Interp qw/:all/;

=head1 DESCRIPTION

Here is a list of all the functions included in this module :

=over 1

=item C<gsl_interp_accel_alloc()> - This function returns a pointer to an accelerator object, which is a kind of iterator for interpolation lookups. It tracks the state of lookups, thus allowing for application of various acceleration strategies.

=item C<gsl_interp_accel_find($a, $x_array, $size, $x)> - This function performs a lookup action on the data array $x_array of size $size, using the given accelerator $a. This is how lookups are performed during evaluation of an interpolation. The function returns an index i such that $x_array[i] <= $x < $x_array[i+1].

=item C<gsl_interp_accel_reset> 

=item C<gsl_interp_accel_free($a)> - This function frees the accelerator object $a. 

=item C<gsl_interp_alloc($T, $alloc)> - This function returns a newly allocated interpolation object of type $T for $size data-points. $T must be one of the constants below.  

=item C<gsl_interp_init($interp, $xa, $ya, $size)> - This function initializes the interpolation object interp for the data (xa,ya) where xa and ya are arrays of size size. The interpolation object (gsl_interp) does not save the data arrays xa and ya and only stores the static state computed from the data. The xa data array is always assumed to be strictly ordered, with increasing x values; the behavior for other arrangements is not defined. 

=item C<gsl_interp_name($interp)> - This function returns the name of the interpolation type used by $interp. 

=item C<gsl_interp_min_size($interp)> - This function returns the minimum number of points required by the interpolation type of $interp. For example, Akima spline interpolation requires a minimum of 5 points.

=item C<gsl_interp_eval_e($interp, $xa, $ya, $x, $acc)> - This functions returns the interpolated value of y for a given point $x, using the interpolation object $interp, data arrays $xa and $ya and the accelerator $acc. The function returns 0 if the operation succeeded, 1 otherwise and the y value.

=item C<gsl_interp_eval($interp, $xa, $ya, $x, $acc)> - This functions returns the interpolated value of y for a given point $x, using the interpolation object $interp, data arrays $xa and $ya and the accelerator $acc. 

=item C<gsl_interp_eval_deriv_e($interp, $xa, $ya, $x, $acc)> - This function computes the derivative value of y for a given point $x, using the interpolation object $interp, data arrays $xa and $ya and the accelerator $acc. The function returns 0 if the operation succeeded, 1 otherwise and the d value.

=item C<gsl_interp_eval_deriv($interp, $xa, $ya, $x, $acc)> - This function returns the derivative d of an interpolated function for a given point $x, using the interpolation object interp, data arrays $xa and $ya and the accelerator $acc.

=item C<gsl_interp_eval_deriv2_e($interp, $xa, $ya, $x, $acc)> - This function computes the second derivative d2 of an interpolated function for a given point $x, using the interpolation object $interp, data arrays $xa and $ya and the accelerator $acc. The function returns 0 if the operation succeeded, 1 otherwise and the d2 value.

=item C<gsl_interp_eval_deriv2($interp, $xa, $ya, $x, $acc)> - This function returns the second derivative d2 of an interpolated function for a given point $x, using the interpolation object $interp, data arrays $xa and $ya and the accelerator $acc.

=item C<gsl_interp_eval_integ_e($interp, $xa, $ya, $a, $b, $acc)> - This function computes the numerical integral result of an interpolated function over the range [$a, $b], using the interpolation object $interp, data arrays $xa and $ya and the accelerator $acc. The function returns 0 if the operation succeeded, 1 otherwise and the result value.

=item C<gsl_interp_eval_integ($interp, $xa, $ya, $a, $b, $acc)> - This function returns the numerical integral result of an interpolated function over the range [$a, $b], using the interpolation object $interp, data arrays $xa and $ya and the accelerator $acc. 

=item C<gsl_interp_free($interp)> - This function frees the interpolation object $interp. 

=item C<gsl_interp_bsearch($x_array, $x, $index_lo, $index_hi)> - This function returns the index i of the array $x_array such that $x_array[i] <= x < $x_array[i+1]. The index is searched for in the range [$index_lo,$index_hi].

=back

This module also includes the following constants :

=over 1

=item C<$gsl_interp_linear> - Linear interpolation

=item C<$gsl_interp_polynomial> - Polynomial interpolation. This method should only be used for interpolating small numbers of points because polynomial interpolation introduces large oscillations, even for well-behaved datasets. The number of terms in the interpolating polynomial is equal to the number of points. 

=item C<$gsl_interp_cspline> - Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second derivatives at the supplied data-points. The second derivative is chosen to be zero at the first point and last point.

=item C<$gsl_interp_cspline_periodic> - Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second derivatives at the supplied data-points. The derivatives at the first and last points are also matched. Note that the last point in the data must have the same y-value as the first point, otherwise the resulting periodic interpolation will have a discontinuity at the boundary. 

=item C<$gsl_interp_akima> - Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded corner algorithm of Wodicka. 

=item C<$gsl_interp_akima_periodic> - Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.

=back

=head1 EXAMPLES

 use Math::GSL::Interp qw/:all/;
 my $x_array = [ 0.0, 1.0, 2.0, 3.0, 4.0 ];
 # check that we get the last interval if x == last value 
 $index_result = gsl_interp_bsearch($x_array, 4.0, 0, 4);
 print "The last interval is $index_result \n";

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut


1;

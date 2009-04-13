# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Deriv;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Derivc;
bootstrap Math::GSL::Deriv;
package Math::GSL::Deriv;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Deriv;

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

package Math::GSL::Deriv;

*gsl_max = *Math::GSL::Derivc::gsl_max;
*gsl_min = *Math::GSL::Derivc::gsl_min;
*gsl_deriv_central = *Math::GSL::Derivc::gsl_deriv_central;
*gsl_deriv_backward = *Math::GSL::Derivc::gsl_deriv_backward;
*gsl_deriv_forward = *Math::GSL::Derivc::gsl_deriv_forward;

############# Class : Math::GSL::Deriv::gsl_function_struct ##############

package Math::GSL::Deriv::gsl_function_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Deriv );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::Derivc::gsl_function_struct_function_get;
*swig_function_set = *Math::GSL::Derivc::gsl_function_struct_function_set;
*swig_params_get = *Math::GSL::Derivc::gsl_function_struct_params_get;
*swig_params_set = *Math::GSL::Derivc::gsl_function_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Derivc::new_gsl_function_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Derivc::delete_gsl_function_struct($self);
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


############# Class : Math::GSL::Deriv::gsl_function_fdf_struct ##############

package Math::GSL::Deriv::gsl_function_fdf_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Deriv );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Derivc::gsl_function_fdf_struct_f_get;
*swig_f_set = *Math::GSL::Derivc::gsl_function_fdf_struct_f_set;
*swig_df_get = *Math::GSL::Derivc::gsl_function_fdf_struct_df_get;
*swig_df_set = *Math::GSL::Derivc::gsl_function_fdf_struct_df_set;
*swig_fdf_get = *Math::GSL::Derivc::gsl_function_fdf_struct_fdf_get;
*swig_fdf_set = *Math::GSL::Derivc::gsl_function_fdf_struct_fdf_set;
*swig_params_get = *Math::GSL::Derivc::gsl_function_fdf_struct_params_get;
*swig_params_set = *Math::GSL::Derivc::gsl_function_fdf_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Derivc::new_gsl_function_fdf_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Derivc::delete_gsl_function_fdf_struct($self);
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


############# Class : Math::GSL::Deriv::gsl_function_vec_struct ##############

package Math::GSL::Deriv::gsl_function_vec_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Deriv );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::Derivc::gsl_function_vec_struct_function_get;
*swig_function_set = *Math::GSL::Derivc::gsl_function_vec_struct_function_set;
*swig_params_get = *Math::GSL::Derivc::gsl_function_vec_struct_params_get;
*swig_params_set = *Math::GSL::Derivc::gsl_function_vec_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Derivc::new_gsl_function_vec_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Derivc::delete_gsl_function_vec_struct($self);
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

package Math::GSL::Deriv;

*GSL_MAJOR_VERSION = *Math::GSL::Derivc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Derivc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Derivc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Derivc::GSL_NEGZERO;
*M_E = *Math::GSL::Derivc::M_E;
*M_LOG2E = *Math::GSL::Derivc::M_LOG2E;
*M_LOG10E = *Math::GSL::Derivc::M_LOG10E;
*M_SQRT2 = *Math::GSL::Derivc::M_SQRT2;
*M_SQRT1_2 = *Math::GSL::Derivc::M_SQRT1_2;
*M_SQRT3 = *Math::GSL::Derivc::M_SQRT3;
*M_PI = *Math::GSL::Derivc::M_PI;
*M_PI_2 = *Math::GSL::Derivc::M_PI_2;
*M_PI_4 = *Math::GSL::Derivc::M_PI_4;
*M_SQRTPI = *Math::GSL::Derivc::M_SQRTPI;
*M_2_SQRTPI = *Math::GSL::Derivc::M_2_SQRTPI;
*M_1_PI = *Math::GSL::Derivc::M_1_PI;
*M_2_PI = *Math::GSL::Derivc::M_2_PI;
*M_LN10 = *Math::GSL::Derivc::M_LN10;
*M_LN2 = *Math::GSL::Derivc::M_LN2;
*M_LNPI = *Math::GSL::Derivc::M_LNPI;
*M_EULER = *Math::GSL::Derivc::M_EULER;

@EXPORT_OK = qw/
               gsl_deriv_central 
               gsl_deriv_backward 
               gsl_deriv_forward 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Deriv - Numerical Derivatives

=head1 SYNOPSIS

    use Math::GSL::Deriv qw/:all/;  
    use Math::GSL::Errno qw/:all/;

    my ($x, $h) = (1.5, 0.01);
    my ($status, $val,$err) = gsl_deriv_central ( sub {  sin($_[0]) }, $x, $h); 
    my $res = abs($val - cos($x));
    if ($status == $GSL_SUCCESS) {
        printf "deriv(sin((%g)) = %.18g, max error=%.18g\n", $x, $val, $err;  
        printf "       cos(%g)) = %.18g, residue=  %.18g\n"  , $x, cos($x), $res;
    } else {
        my $gsl_error = gsl_strerror($status);
        print "Numerical Derivative FAILED, reason:\n $gsl_error\n\n";
    }


=head1 DESCRIPTION

This module allows you to take the numerical derivative of a Perl subroutine. To find
a numerical derivative you must also specify a point to evaluate the derivative and a
"step size". The step size is a knob that you can turn to get a more finely or coarse
grained approximation. As the step size $h goes to zero, the formal definition of a
derivative is reached, but in practive you must choose a reasonable step size to get
a reasonable answer. Usually something in the range of 1/10 to 1/10000 is sufficient.

So long as your function returns a single scalar value, you can differentiate as 
complicated a function as your heart desires.

=over 

=item * C<gsl_deriv_central($function, $x, $h)> 

    use Math::GSL::Deriv qw/gsl_deriv_central/;
    my ($x, $h) = (1.5, 0.01);
    sub func { my $x=shift; $x**4 - 15 * $x + sqrt($x) };

    my ($status, $val,$err) = gsl_deriv_central ( \&func , $x, $h); 

    This method approximates the central difference of the subroutine reference
    $function, evaluated at $x, with "step size" $h. This means that the
    function is evaluated at $x-$h and $x+h.


=item * C<gsl_deriv_backward($function, $x, $h)> 

    use Math::GSL::Deriv qw/gsl_deriv_backward/;
    my ($x, $h) = (1.5, 0.01);
    sub func { my $x=shift; $x**4 - 15 * $x + sqrt($x) };

    my ($status, $val,$err) = gsl_deriv_backward ( \&func , $x, $h); 

    This method approximates the backward difference of the subroutine
    reference $function, evaluated at $x, with "step size" $h. This means that
    the function is evaluated at $x-$h and $x.

=item * C<gsl_deriv_forward($function, $x, $h)> 

    use Math::GSL::Deriv qw/gsl_deriv_forward/;
    my ($x, $h) = (1.5, 0.01);
    sub func { my $x=shift; $x**4 - 15 * $x + sqrt($x) };

    my ($status, $val,$err) = gsl_deriv_forward ( \&func , $x, $h); 

    This method approximates the forward difference of the subroutine reference
    $function, evaluated at $x, with "step size" $h. This means that the
    function is evaluated at $x and $x+$h.

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

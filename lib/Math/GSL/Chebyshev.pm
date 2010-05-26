# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.40
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Chebyshev;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Chebyshevc;
bootstrap Math::GSL::Chebyshev;
package Math::GSL::Chebyshev;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Chebyshev;

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

package Math::GSL::Chebyshev;

*gsl_cheb_alloc = *Math::GSL::Chebyshevc::gsl_cheb_alloc;
*gsl_cheb_free = *Math::GSL::Chebyshevc::gsl_cheb_free;
*gsl_cheb_init = *Math::GSL::Chebyshevc::gsl_cheb_init;
*gsl_cheb_order = *Math::GSL::Chebyshevc::gsl_cheb_order;
*gsl_cheb_size = *Math::GSL::Chebyshevc::gsl_cheb_size;
*gsl_cheb_coeffs = *Math::GSL::Chebyshevc::gsl_cheb_coeffs;
*gsl_cheb_eval = *Math::GSL::Chebyshevc::gsl_cheb_eval;
*gsl_cheb_eval_err = *Math::GSL::Chebyshevc::gsl_cheb_eval_err;
*gsl_cheb_eval_n = *Math::GSL::Chebyshevc::gsl_cheb_eval_n;
*gsl_cheb_eval_n_err = *Math::GSL::Chebyshevc::gsl_cheb_eval_n_err;
*gsl_cheb_eval_mode = *Math::GSL::Chebyshevc::gsl_cheb_eval_mode;
*gsl_cheb_eval_mode_e = *Math::GSL::Chebyshevc::gsl_cheb_eval_mode_e;
*gsl_cheb_calc_deriv = *Math::GSL::Chebyshevc::gsl_cheb_calc_deriv;
*gsl_cheb_calc_integ = *Math::GSL::Chebyshevc::gsl_cheb_calc_integ;

############# Class : Math::GSL::Chebyshev::gsl_cheb_series_struct ##############

package Math::GSL::Chebyshev::gsl_cheb_series_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Chebyshev );
%OWNER = ();
%ITERATORS = ();
*swig_c_get = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_c_get;
*swig_c_set = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_c_set;
*swig_order_get = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_order_get;
*swig_order_set = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_order_set;
*swig_a_get = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_a_get;
*swig_a_set = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_a_set;
*swig_b_get = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_b_get;
*swig_b_set = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_b_set;
*swig_order_sp_get = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_order_sp_get;
*swig_order_sp_set = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_order_sp_set;
*swig_f_get = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_f_get;
*swig_f_set = *Math::GSL::Chebyshevc::gsl_cheb_series_struct_f_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Chebyshevc::new_gsl_cheb_series_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Chebyshevc::delete_gsl_cheb_series_struct($self);
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


############# Class : Math::GSL::Chebyshev::gsl_function_struct ##############

package Math::GSL::Chebyshev::gsl_function_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Chebyshev );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::Chebyshevc::gsl_function_struct_function_get;
*swig_function_set = *Math::GSL::Chebyshevc::gsl_function_struct_function_set;
*swig_params_get = *Math::GSL::Chebyshevc::gsl_function_struct_params_get;
*swig_params_set = *Math::GSL::Chebyshevc::gsl_function_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Chebyshevc::new_gsl_function_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Chebyshevc::delete_gsl_function_struct($self);
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


############# Class : Math::GSL::Chebyshev::gsl_function_fdf_struct ##############

package Math::GSL::Chebyshev::gsl_function_fdf_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Chebyshev );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_f_get;
*swig_f_set = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_f_set;
*swig_df_get = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_df_get;
*swig_df_set = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_df_set;
*swig_fdf_get = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_fdf_get;
*swig_fdf_set = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_fdf_set;
*swig_params_get = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_params_get;
*swig_params_set = *Math::GSL::Chebyshevc::gsl_function_fdf_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Chebyshevc::new_gsl_function_fdf_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Chebyshevc::delete_gsl_function_fdf_struct($self);
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


############# Class : Math::GSL::Chebyshev::gsl_function_vec_struct ##############

package Math::GSL::Chebyshev::gsl_function_vec_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Chebyshev );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::Chebyshevc::gsl_function_vec_struct_function_get;
*swig_function_set = *Math::GSL::Chebyshevc::gsl_function_vec_struct_function_set;
*swig_params_get = *Math::GSL::Chebyshevc::gsl_function_vec_struct_params_get;
*swig_params_set = *Math::GSL::Chebyshevc::gsl_function_vec_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Chebyshevc::new_gsl_function_vec_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Chebyshevc::delete_gsl_function_vec_struct($self);
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

package Math::GSL::Chebyshev;

*GSL_MAJOR_VERSION = *Math::GSL::Chebyshevc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Chebyshevc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Chebyshevc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Chebyshevc::GSL_NEGZERO;
*M_E = *Math::GSL::Chebyshevc::M_E;
*M_LOG2E = *Math::GSL::Chebyshevc::M_LOG2E;
*M_LOG10E = *Math::GSL::Chebyshevc::M_LOG10E;
*M_SQRT2 = *Math::GSL::Chebyshevc::M_SQRT2;
*M_SQRT1_2 = *Math::GSL::Chebyshevc::M_SQRT1_2;
*M_SQRT3 = *Math::GSL::Chebyshevc::M_SQRT3;
*M_PI = *Math::GSL::Chebyshevc::M_PI;
*M_PI_2 = *Math::GSL::Chebyshevc::M_PI_2;
*M_PI_4 = *Math::GSL::Chebyshevc::M_PI_4;
*M_SQRTPI = *Math::GSL::Chebyshevc::M_SQRTPI;
*M_2_SQRTPI = *Math::GSL::Chebyshevc::M_2_SQRTPI;
*M_1_PI = *Math::GSL::Chebyshevc::M_1_PI;
*M_2_PI = *Math::GSL::Chebyshevc::M_2_PI;
*M_LN10 = *Math::GSL::Chebyshevc::M_LN10;
*M_LN2 = *Math::GSL::Chebyshevc::M_LN2;
*M_LNPI = *Math::GSL::Chebyshevc::M_LNPI;
*M_EULER = *Math::GSL::Chebyshevc::M_EULER;
*GSL_PREC_DOUBLE = *Math::GSL::Chebyshevc::GSL_PREC_DOUBLE;
*GSL_PREC_SINGLE = *Math::GSL::Chebyshevc::GSL_PREC_SINGLE;
*GSL_PREC_APPROX = *Math::GSL::Chebyshevc::GSL_PREC_APPROX;
*GSL_MODE_DEFAULT = *Math::GSL::Chebyshevc::GSL_MODE_DEFAULT;

@EXPORT_OK = qw/
               gsl_cheb_alloc 
               gsl_cheb_free 
               gsl_cheb_init 
               gsl_cheb_eval 
               gsl_cheb_eval_err 
               gsl_cheb_eval_n 
               gsl_cheb_eval_n_err 
               gsl_cheb_eval_mode 
               gsl_cheb_eval_mode_e 
               gsl_cheb_calc_deriv 
               gsl_cheb_calc_integ 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Chebyshev - Univariate Chebyshev Series Approximation

=head1 SYNOPSIS

    use Math::GSL::Chebyshev qw /:all/;

    my $cheb             = gsl_cheb_alloc(40);
    my $function         = sub { sin(cos($_[0])) };

    gsl_cheb_init($cheb, $function, 0, 10);

    my $x                = gsl_cheb_eval($cheb, 5.5 );
    my ($status,$y,$err) = gsl_cheb_eval_err($cheb, 7.5 );
    gsl_cheb_free($cheb);

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over

=item * C<gsl_cheb_alloc($size)>

    my $cheb = gsl_cheb_alloc(50);

Allocates a new Chebyshev object with $size sample points.

=item * C<gsl_cheb_free($cheb)>

Deallocates memory associated to $cheb. Returns void.

=item * C<gsl_cheb_init($cheb,$function, $lower, $upper)>

    gsl_cheb_init($cheb, sub { sin(cos($_[0])) }, 0, 10 );

Initiate a Chebyshev object with a function and upper and lower bounds.
Returns void.

=item * C<gsl_cheb_eval($function, $value)>

    my $evaluated = gsl_cheb_eval($cheb, 5 );

Returns a Perl scalar of the Chebyshev object $cheb evaluated at $value.

=item * C<gsl_cheb_eval_err($cheb, $value)>

    my ($status,$evaluated,$err) = gsl_cheb_eval($cheb, 5 );

Returns a list consisting of a GSL status code, the evaluate value and
the estimated error of the evaluation.

=item * C<gsl_cheb_eval_n >

=item * C<gsl_cheb_eval_n_err >

=item * C<gsl_cheb_eval_mode >

=item * C<gsl_cheb_eval_mode_e >

=item * C<gsl_cheb_calc_deriv($deriv,$cheb) >

   my $status = gsl_cheb_calc_deriv( $deriv, $cheb ); 

This will calculate the derivative of $cheb and stores it
in $deriv, which must be pre-allocated. Returns a GSL status code.

=item * C<gsl_cheb_calc_integ($integ,$cheb) >

   my $status = gsl_cheb_calc_integ( $deriv, $cheb ); 

This will calculate the derivative of $cheb and stores it
in $deriv, which must be pre-allocated. Returns a GSL status code.

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

 Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want


=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

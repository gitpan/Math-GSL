# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Complex;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Complexc;
bootstrap Math::GSL::Complex;
package Math::GSL::Complex;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Complex;

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

package Math::GSL::Complex;

*gsl_complex_polar = *Math::GSL::Complexc::gsl_complex_polar;
*gsl_complex_rect = *Math::GSL::Complexc::gsl_complex_rect;
*gsl_complex_arg = *Math::GSL::Complexc::gsl_complex_arg;
*gsl_complex_abs = *Math::GSL::Complexc::gsl_complex_abs;
*gsl_complex_abs2 = *Math::GSL::Complexc::gsl_complex_abs2;
*gsl_complex_logabs = *Math::GSL::Complexc::gsl_complex_logabs;
*gsl_complex_add = *Math::GSL::Complexc::gsl_complex_add;
*gsl_complex_sub = *Math::GSL::Complexc::gsl_complex_sub;
*gsl_complex_mul = *Math::GSL::Complexc::gsl_complex_mul;
*gsl_complex_div = *Math::GSL::Complexc::gsl_complex_div;
*gsl_complex_add_real = *Math::GSL::Complexc::gsl_complex_add_real;
*gsl_complex_sub_real = *Math::GSL::Complexc::gsl_complex_sub_real;
*gsl_complex_mul_real = *Math::GSL::Complexc::gsl_complex_mul_real;
*gsl_complex_div_real = *Math::GSL::Complexc::gsl_complex_div_real;
*gsl_complex_add_imag = *Math::GSL::Complexc::gsl_complex_add_imag;
*gsl_complex_sub_imag = *Math::GSL::Complexc::gsl_complex_sub_imag;
*gsl_complex_mul_imag = *Math::GSL::Complexc::gsl_complex_mul_imag;
*gsl_complex_div_imag = *Math::GSL::Complexc::gsl_complex_div_imag;
*gsl_complex_conjugate = *Math::GSL::Complexc::gsl_complex_conjugate;
*gsl_complex_inverse = *Math::GSL::Complexc::gsl_complex_inverse;
*gsl_complex_negative = *Math::GSL::Complexc::gsl_complex_negative;
*gsl_complex_sqrt = *Math::GSL::Complexc::gsl_complex_sqrt;
*gsl_complex_sqrt_real = *Math::GSL::Complexc::gsl_complex_sqrt_real;
*gsl_complex_pow = *Math::GSL::Complexc::gsl_complex_pow;
*gsl_complex_pow_real = *Math::GSL::Complexc::gsl_complex_pow_real;
*gsl_complex_exp = *Math::GSL::Complexc::gsl_complex_exp;
*gsl_complex_log = *Math::GSL::Complexc::gsl_complex_log;
*gsl_complex_log10 = *Math::GSL::Complexc::gsl_complex_log10;
*gsl_complex_log_b = *Math::GSL::Complexc::gsl_complex_log_b;
*gsl_complex_sin = *Math::GSL::Complexc::gsl_complex_sin;
*gsl_complex_cos = *Math::GSL::Complexc::gsl_complex_cos;
*gsl_complex_sec = *Math::GSL::Complexc::gsl_complex_sec;
*gsl_complex_csc = *Math::GSL::Complexc::gsl_complex_csc;
*gsl_complex_tan = *Math::GSL::Complexc::gsl_complex_tan;
*gsl_complex_cot = *Math::GSL::Complexc::gsl_complex_cot;
*gsl_complex_arcsin = *Math::GSL::Complexc::gsl_complex_arcsin;
*gsl_complex_arcsin_real = *Math::GSL::Complexc::gsl_complex_arcsin_real;
*gsl_complex_arccos = *Math::GSL::Complexc::gsl_complex_arccos;
*gsl_complex_arccos_real = *Math::GSL::Complexc::gsl_complex_arccos_real;
*gsl_complex_arcsec = *Math::GSL::Complexc::gsl_complex_arcsec;
*gsl_complex_arcsec_real = *Math::GSL::Complexc::gsl_complex_arcsec_real;
*gsl_complex_arccsc = *Math::GSL::Complexc::gsl_complex_arccsc;
*gsl_complex_arccsc_real = *Math::GSL::Complexc::gsl_complex_arccsc_real;
*gsl_complex_arctan = *Math::GSL::Complexc::gsl_complex_arctan;
*gsl_complex_arccot = *Math::GSL::Complexc::gsl_complex_arccot;
*gsl_complex_sinh = *Math::GSL::Complexc::gsl_complex_sinh;
*gsl_complex_cosh = *Math::GSL::Complexc::gsl_complex_cosh;
*gsl_complex_sech = *Math::GSL::Complexc::gsl_complex_sech;
*gsl_complex_csch = *Math::GSL::Complexc::gsl_complex_csch;
*gsl_complex_tanh = *Math::GSL::Complexc::gsl_complex_tanh;
*gsl_complex_coth = *Math::GSL::Complexc::gsl_complex_coth;
*gsl_complex_arcsinh = *Math::GSL::Complexc::gsl_complex_arcsinh;
*gsl_complex_arccosh = *Math::GSL::Complexc::gsl_complex_arccosh;
*gsl_complex_arccosh_real = *Math::GSL::Complexc::gsl_complex_arccosh_real;
*gsl_complex_arcsech = *Math::GSL::Complexc::gsl_complex_arcsech;
*gsl_complex_arccsch = *Math::GSL::Complexc::gsl_complex_arccsch;
*gsl_complex_arctanh = *Math::GSL::Complexc::gsl_complex_arctanh;
*gsl_complex_arctanh_real = *Math::GSL::Complexc::gsl_complex_arctanh_real;
*gsl_complex_arccoth = *Math::GSL::Complexc::gsl_complex_arccoth;
*new_doubleArray = *Math::GSL::Complexc::new_doubleArray;
*delete_doubleArray = *Math::GSL::Complexc::delete_doubleArray;
*doubleArray_getitem = *Math::GSL::Complexc::doubleArray_getitem;
*doubleArray_setitem = *Math::GSL::Complexc::doubleArray_setitem;

############# Class : Math::GSL::Complex::gsl_complex_long_double ##############

package Math::GSL::Complex::gsl_complex_long_double;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Complex );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::Complexc::gsl_complex_long_double_dat_get;
*swig_dat_set = *Math::GSL::Complexc::gsl_complex_long_double_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Complexc::new_gsl_complex_long_double(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Complexc::delete_gsl_complex_long_double($self);
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


############# Class : Math::GSL::Complex::gsl_complex ##############

package Math::GSL::Complex::gsl_complex;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Complex );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::Complexc::gsl_complex_dat_get;
*swig_dat_set = *Math::GSL::Complexc::gsl_complex_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Complexc::new_gsl_complex(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Complexc::delete_gsl_complex($self);
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


############# Class : Math::GSL::Complex::gsl_complex_float ##############

package Math::GSL::Complex::gsl_complex_float;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Complex );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::Complexc::gsl_complex_float_dat_get;
*swig_dat_set = *Math::GSL::Complexc::gsl_complex_float_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Complexc::new_gsl_complex_float(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Complexc::delete_gsl_complex_float($self);
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

package Math::GSL::Complex;

*GSL_MAJOR_VERSION = *Math::GSL::Complexc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Complexc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Complexc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Complexc::GSL_NEGZERO;


@EXPORT_OK = qw(
    gsl_complex_arg gsl_complex_abs gsl_complex_rect gsl_complex_polar doubleArray_getitem 
    gsl_complex_rect gsl_complex_polar gsl_complex_arg gsl_complex_abs gsl_complex_abs2 
    gsl_complex_logabs gsl_complex_add gsl_complex_sub gsl_complex_mul gsl_complex_div 
    gsl_complex_add_real gsl_complex_sub_real gsl_complex_mul_real gsl_complex_div_real 
    gsl_complex_add_imag gsl_complex_sub_imag gsl_complex_mul_imag gsl_complex_div_imag 
    gsl_complex_conjugate gsl_complex_inverse gsl_complex_negative gsl_complex_sqrt 
    gsl_complex_sqrt_real gsl_complex_pow gsl_complex_pow_real gsl_complex_exp 
    gsl_complex_log gsl_complex_log10 gsl_complex_log_b gsl_complex_sin 
    gsl_complex_cos gsl_complex_sec gsl_complex_csc gsl_complex_tan 
    gsl_complex_cot gsl_complex_arcsin gsl_complex_arcsin_real gsl_complex_arccos 
    gsl_complex_arccos_real gsl_complex_arcsec gsl_complex_arcsec_real gsl_complex_arccsc 
    gsl_complex_arccsc_real gsl_complex_arctan gsl_complex_arccot gsl_complex_sinh 
    gsl_complex_cosh gsl_complex_sech gsl_complex_csch gsl_complex_tanh 
    gsl_complex_coth gsl_complex_arcsinh gsl_complex_arccosh gsl_complex_arccosh_real 
    gsl_complex_arcsech gsl_complex_arccsch gsl_complex_arctanh gsl_complex_arctanh_real 
    gsl_complex_arccoth new_doubleArray delete_doubleArray doubleArray_setitem
    gsl_real gsl_imag gsl_parts
    gsl_complex_eq gsl_set_real gsl_set_imag gsl_set_complex
    $GSL_COMPLEX_ONE $GSL_COMPLEX_ZERO $GSL_COMPLEX_NEGONE
);
# macros to implement
# gsl_set_complex gsl_set_complex_packed
our ($GSL_COMPLEX_ONE, $GSL_COMPLEX_ZERO, $GSL_COMPLEX_NEGONE) = map { gsl_complex_rect($_, 0) } qw(1 0 -1); 


%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

sub new {
    my ($class, @values) = @_;
    my $this = {};
    $this->{_complex} = gsl_complex_rect($values[0], $values[1]);
    bless $this, $class;
}
sub real {
    my ($self) = @_;
    gsl_real($self->{_complex}->{dat});
}

sub imag {
    my ($self) = @_;
    gsl_imag($self->{_complex}->{dat});
}

sub parts {
    my ($self) = @_;
    gsl_parts($self->{_complex}->{dat});
}

sub raw  { (shift)->{_complex} }


### some important macros that are in gsl_complex.h
sub gsl_complex_eq {
    my ($z,$w) = @_;
    gsl_real($z) == gsl_real($w) && gsl_imag($z) == gsl_imag($w) ? 1 : 0;
}

sub gsl_set_real {
    my ($z,$r) = @_;
    doubleArray_setitem($z->{dat}, 0, $r);
}

sub gsl_set_imag {
    my ($z,$i) = @_;
    doubleArray_setitem($z->{dat}, 1, $i);
}

sub gsl_real {
    my $z = shift;
    return doubleArray_getitem($z->{dat}, 0 );
}

sub gsl_imag {
    my $z = shift;
    return doubleArray_getitem($z->{dat}, 1 );
}

sub gsl_parts {
    my $z = shift;
    return (gsl_real($z), gsl_imag($z));
}

sub gsl_set_complex {
    my ($z, $r, $i) = @_;
    gsl_set_real($z, $r);
    gsl_set_imag($z, $i);
}

=head1 NAME

Math::GSL::Complex - Complex Numbers

=head1 SYNOPSIS

    use Math::GSL::Complex qw/:all/;
    my $complex = Math::GSL::Complex->new([3,2]); # creates a complex number 3+2*i
    my $real = $complex->real;                    # returns the real part
    my $imag = $complex->imag;                    # returns the imaginary part
    $complex->gsl_set_real(5);                    # changes the real part to 5
    $complex->gsl_set_imag(4);                    # changes the imaginary part to 4
    $complex->gsl_set_complex(7,6);               # changes it to 7 + 6*I
    ($real, $imag) = $complex->parts;

=head1 DESCRIPTION

Here is a list of all the functions included in this module :

=over 1

=item C<gsl_complex_arg($z)> - return the argument of the complex number $z 

=item C<gsl_complex_abs($z)> - return |$z|, the magnitude of the complex number $z 

=item C<gsl_complex_rect($x,$y)> - create a complex number in cartesian form $x + $y*I

=item C<gsl_complex_polar($r,$theta)> - create a complex number in polar form $r*exp(I*$theta) 

=item C<gsl_complex_abs2($z)> - return |$z|^2, the squared magnitude of the complex number $z

=item C<gsl_complex_logabs($z)> - return log(|$z|), the natural logarithm of the magnitude of the complex number $z 

=item C<gsl_complex_add($c1, $c2)> - return a complex number which is the sum of the complex numbers $c1 and $c2 

=item C<gsl_complex_sub($c1, $c2)> - return a complex number which is the difference between $c1 and $c2 ($c1 - $c2) 

=item C<gsl_complex_mul($c1, $c2)> - return a complex number which is the product of the complex numbers $c1 and $c2

=item C<gsl_complex_div($c1, $c2)> - return a complex number which is the quotient of the complex numbers $c1 and $c2 ($c1 / $c2)

=item C<gsl_complex_add_real($c, $x)> - return the sum of the complex number $c and the real number $x 

=item C<gsl_complex_sub_real($c, $x)> - return the difference of the complex number $c and the real number $x 

=item C<gsl_complex_mul_real($c, $x)> - return the product of the complex number $c and the real number $x 

=item C<gsl_complex_div_real($c, $x)> - return the quotient of the complex number $c and the real number $x  

=item C<gsl_complex_add_imag($c, $y)> - return sum of the complex number $c and the imaginary number i*$x 

=item C<gsl_complex_sub_imag($c, $y)> - return the diffrence of the complex number $c and the imaginary number i*$x 

=item C<gsl_complex_mul_imag($c, $y)> - return the product of the complex number $c and the imaginary number i*$x  

=item C<gsl_complex_div_imag($c, $y)> - return the quotient of the complex number $c and the imaginary number i*$x 

=item C<gsl_complex_conjugate($c)> - return the conjugate of the of the complex number $c (x - i*y)  

=item C<gsl_complex_inverse($c)> - return the inverse, or reciprocal of the complex number $c (1/$c) 

=item C<gsl_complex_negative($c)> - return the negative of the complex number $c (-x -i*y) 

=item C<gsl_complex_sqrt($c)> - return the square root of the complex number $c 

=item C<gsl_complex_sqrt_real($x)> - return the complex square root of the real number $x, where $x may be negative

=item C<gsl_complex_pow($c1, $c2)> - return the complex number $c1 raised to the complex power $c2 

=item C<gsl_complex_pow_real($c, $x)> - return the complex number raised to the real power $x 

=item C<gsl_complex_exp($c)> - return the complex exponential of the complex number $c 

=item C<gsl_complex_log($c)> - return the complex natural logarithm (base e) of the complex number $c 

=item C<gsl_complex_log10($c)> - return the complex base-10 logarithm of the complex number $c

=item C<gsl_complex_log_b($c, $b)> - return the complex base-$b of the complex number $c 

=item C<gsl_complex_sin($c)> - return the complex sine of the complex number $c

=item C<gsl_complex_cos($c)> - return the complex cosine of the complex number $c 

=item C<gsl_complex_sec($c)> - return the complex secant of the complex number $c 

=item C<gsl_complex_csc($c)> - return the complex cosecant of the complex number $c 

=item C<gsl_complex_tan($c)> - return the complex tangent of the complex number $c 

=item C<gsl_complex_cot($c)> - return the complex cotangent of the complex number $c 

=item C<gsl_complex_arcsin($c)> - return the complex arcsine of the complex number $c 

=item C<gsl_complex_arcsin_real($x)> - return the complex arcsine of the real number $x 

=item C<gsl_complex_arccos($c)> - return the complex arccosine of the complex number $c 

=item C<gsl_complex_arccos_real($x)> - return the complex arccosine of the real number $x 

=item C<gsl_complex_arcsec($c)> - return the complex arcsecant of the complex number $c 

=item C<gsl_complex_arcsec_real($x)> - return the complex arcsecant of the real number $x

=item C<gsl_complex_arccsc($c)> - return the complex arccosecant of the complex number $c 

=item C<gsl_complex_arccsc_real($x)> - return the complex arccosecant of the real number $x

=item C<gsl_complex_arctan($c)> - return the complex arctangent of the complex number $c

=item C<gsl_complex_arccot($c)> - return the complex arccotangent of the complex number $c 

=item C<gsl_complex_sinh($c)> - return the complex hyperbolic sine of the complex number $c 

=item C<gsl_complex_cosh($c)> - return the complex hyperbolic cosine of the complex number $cy

=item C<gsl_complex_sech($c)> - return the complex hyperbolic secant of the complex number $c

=item C<gsl_complex_csch($c)> - return the complex hyperbolic cosecant of the complex number $c

=item C<gsl_complex_tanh($c)> - return the complex hyperbolic tangent of the complex number $c

=item C<gsl_complex_coth($c)> - return the complex hyperbolic cotangent of the complex number $c

=item C<gsl_complex_arcsinh($c)> - return the complex hyperbolic arcsine of the complex number $c

=item C<gsl_complex_arccosh($c)> - return the complex hyperbolic arccosine of the complex number $c

=item C<gsl_complex_arccosh_real($x)> - return the complex hyperbolic arccosine of the real number $x 

=item C<gsl_complex_arcsech($c)> - return the complex hyperbolic arcsecant of the complex number $c

=item C<gsl_complex_arccsch($c)> - return the complex hyperbolic arccosecant of the complex number $c

=item C<gsl_complex_arctanh($c)> - return the complex hyperbolic arctangent of the complex number $c

=item C<gsl_complex_arctanh_real($x)> - return the complex hyperbolic arctangent of the real number $x 

=item C<gsl_complex_arccoth($c)> - return the complex hyperbolic arccotangent of the complex number $c

=item C<gsl_real($z)> - return the real part of $z 

=item C<gsl_imag($z)> - return the imaginary part of $z 

=item C<gsl_parts($z)> - return a list of the real and imaginary parts of $z

=item C<gsl_set_real($z, $x)> - sets the real part of $z to $x

=item C<gsl_set_imag($z, $y)> - sets the imaginary part of $z to $y

=item C<gsl_set_complex($z, $x, $h)> - sets the real part of $z to $x and the imaginary part to $y

=back

=head1 EXAMPLES

This code defines $z as 6 + 4*I, takes the complex conjugate of that number, then prints it out.

=over 1

    my $z = gsl_complex_rect(6,4);
    my $y = gsl_complex_conjugate($z);
    my ($real, $imag) = gsl_parts($y);
    print "z = $real + $imag*I\n";

=back

This code defines $z as 5 + 3*I, multiplies it by 2 and then prints it out.

=over 1

    my $x = gsl_complex_rect(5,3);
    my $z = gsl_complex_mul_real($x, 2);
    my $real = gsl_real($z);
    my $imag = gsl_imag($z);
    print "Re(\$z) = $real\n";

=back

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
1;

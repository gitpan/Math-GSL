# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::PowInt;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::PowIntc;
bootstrap Math::GSL::PowInt;
package Math::GSL::PowInt;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::PowInt;

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

package Math::GSL::PowInt;

*gsl_pow_2 = *Math::GSL::PowIntc::gsl_pow_2;
*gsl_pow_3 = *Math::GSL::PowIntc::gsl_pow_3;
*gsl_pow_4 = *Math::GSL::PowIntc::gsl_pow_4;
*gsl_pow_5 = *Math::GSL::PowIntc::gsl_pow_5;
*gsl_pow_6 = *Math::GSL::PowIntc::gsl_pow_6;
*gsl_pow_7 = *Math::GSL::PowIntc::gsl_pow_7;
*gsl_pow_8 = *Math::GSL::PowIntc::gsl_pow_8;
*gsl_pow_9 = *Math::GSL::PowIntc::gsl_pow_9;
*gsl_pow_int = *Math::GSL::PowIntc::gsl_pow_int;

# ------- VARIABLE STUBS --------

package Math::GSL::PowInt;

*GSL_MAJOR_VERSION = *Math::GSL::PowIntc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::PowIntc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::PowIntc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::PowIntc::GSL_NEGZERO;


our @EXPORT_OK = qw/gsl_pow_2 gsl_pow_2 gsl_pow_3 
                    gsl_pow_4 gsl_pow_5 gsl_pow_6 
                    gsl_pow_7 gsl_pow_8 gsl_pow_9 gsl_pow_int
                /;
our %EXPORT_TAGS = ( all =>  \@EXPORT_OK );

__END__

=head1 NAME

Math::GSL::PowInt - Integer Power functions

=head1 SYNOPSIS

    use Math::GSL::PowInt qw /gsl_pow_2 gsl_pow_4 gsl_pow_int/;
    print '2**4  = ' . gsl_pow_2(4) . "\n";
    print '4**7  = ' . gsl_pow_4(7) . "\n";
    print '17**5 = ' . gsl_pow_int(17,5) . "\n";

=head1 DESCRIPTION

This module implements the GSL Integer Power functions, which allow one to
quickly compute a given integer raised to any real number.  It contains
gsl_pow_2 to gsl_pow_9 and gsl_pow_int. If you need a power higher than 9, you
can use gsl_pow_int, which takes a base as the first argument and power being
raised to as the second argument.

You can also write 

C<use Math::GSL::PowInt qw/:all/;>
    
to use all avaible functions of the module.

=head1 PURPOSE

This module is used to speed up arithmetic operations. The underlying code
decomposes the multiplication in the most efficient way. Since it is
implemented in XS (via SWIG), it should lend a large performance increase over
the Perl builtin exponentiation operator, '**' .

=head1 BENCHMARKS

Would someone like to submit some benchmarks?

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
1;

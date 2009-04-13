# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Spline;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Splinec;
bootstrap Math::GSL::Spline;
package Math::GSL::Spline;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Spline;

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

package Math::GSL::Spline;

*gsl_spline_alloc = *Math::GSL::Splinec::gsl_spline_alloc;
*gsl_spline_init = *Math::GSL::Splinec::gsl_spline_init;
*gsl_spline_name = *Math::GSL::Splinec::gsl_spline_name;
*gsl_spline_min_size = *Math::GSL::Splinec::gsl_spline_min_size;
*gsl_spline_eval_e = *Math::GSL::Splinec::gsl_spline_eval_e;
*gsl_spline_eval = *Math::GSL::Splinec::gsl_spline_eval;
*gsl_spline_eval_deriv_e = *Math::GSL::Splinec::gsl_spline_eval_deriv_e;
*gsl_spline_eval_deriv = *Math::GSL::Splinec::gsl_spline_eval_deriv;
*gsl_spline_eval_deriv2_e = *Math::GSL::Splinec::gsl_spline_eval_deriv2_e;
*gsl_spline_eval_deriv2 = *Math::GSL::Splinec::gsl_spline_eval_deriv2;
*gsl_spline_eval_integ_e = *Math::GSL::Splinec::gsl_spline_eval_integ_e;
*gsl_spline_eval_integ = *Math::GSL::Splinec::gsl_spline_eval_integ;
*gsl_spline_free = *Math::GSL::Splinec::gsl_spline_free;

############# Class : Math::GSL::Spline::gsl_spline ##############

package Math::GSL::Spline::gsl_spline;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Spline );
%OWNER = ();
%ITERATORS = ();
*swig_interp_get = *Math::GSL::Splinec::gsl_spline_interp_get;
*swig_interp_set = *Math::GSL::Splinec::gsl_spline_interp_set;
*swig_x_get = *Math::GSL::Splinec::gsl_spline_x_get;
*swig_x_set = *Math::GSL::Splinec::gsl_spline_x_set;
*swig_y_get = *Math::GSL::Splinec::gsl_spline_y_get;
*swig_y_set = *Math::GSL::Splinec::gsl_spline_y_set;
*swig_size_get = *Math::GSL::Splinec::gsl_spline_size_get;
*swig_size_set = *Math::GSL::Splinec::gsl_spline_size_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Splinec::new_gsl_spline(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Splinec::delete_gsl_spline($self);
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

package Math::GSL::Spline;

*GSL_MAJOR_VERSION = *Math::GSL::Splinec::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Splinec::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Splinec::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Splinec::GSL_NEGZERO;

@EXPORT_OK = qw/
               gsl_spline_alloc 
               gsl_spline_init 
               gsl_spline_name 
               gsl_spline_min_size 
               gsl_spline_eval_e 
               gsl_spline_eval 
               gsl_spline_eval_deriv_e 
               gsl_spline_eval_deriv 
               gsl_spline_eval_deriv2_e 
               gsl_spline_eval_deriv2 
               gsl_spline_eval_integ_e 
               gsl_spline_eval_integ 
               gsl_spline_free 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Spline - Higher-level Interface to Interp

=head1 SYNOPSIS

use Math::GSL::Spline qw /:all/;

=head1 DESCRIPTION

The functions described in the Interp module required the user to supply pointers to the x and y arrays on each call. The following functions are equivalent to the corresponding gsl_interp functions but maintain a copy of this data in the gsl_spline object. This removes the need to pass both xa and ya as arguments on each evaluation.

Here is a list of all the functions in this module :

=over

=item * C<gsl_spline_alloc($T, $size)>

=item * C<gsl_spline_init($spline, $xa, $ya, $size)>

=item * C<gsl_spline_free($spline)>

=item * C<gsl_spline_name($spline)>

=item * C<gsl_spline_min_size($spline)>

=item * C<gsl_spline_eval_e($spline, $x, $acc)>

=item * C<gsl_spline_eval($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv_e($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv2_e($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv2($spline, $x, $acc)>

=item * C<gsl_spline_eval_integ_e($spline, $a, $b, $acc)>

=item * C<gsl_spline_eval_integ($spline, $a, $b, $acc)>

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

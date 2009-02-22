# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::QRNG;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::QRNGc;
bootstrap Math::GSL::QRNG;
package Math::GSL::QRNG;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::QRNG;

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

package Math::GSL::QRNG;

*gsl_qrng_alloc = *Math::GSL::QRNGc::gsl_qrng_alloc;
*gsl_qrng_memcpy = *Math::GSL::QRNGc::gsl_qrng_memcpy;
*gsl_qrng_clone = *Math::GSL::QRNGc::gsl_qrng_clone;
*gsl_qrng_free = *Math::GSL::QRNGc::gsl_qrng_free;
*gsl_qrng_init = *Math::GSL::QRNGc::gsl_qrng_init;
*gsl_qrng_name = *Math::GSL::QRNGc::gsl_qrng_name;
*gsl_qrng_size = *Math::GSL::QRNGc::gsl_qrng_size;
*gsl_qrng_state = *Math::GSL::QRNGc::gsl_qrng_state;
*gsl_qrng_get = *Math::GSL::QRNGc::gsl_qrng_get;

############# Class : Math::GSL::QRNG::gsl_qrng_type ##############

package Math::GSL::QRNG::gsl_qrng_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::QRNG );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::QRNGc::gsl_qrng_type_name_get;
*swig_name_set = *Math::GSL::QRNGc::gsl_qrng_type_name_set;
*swig_max_dimension_get = *Math::GSL::QRNGc::gsl_qrng_type_max_dimension_get;
*swig_max_dimension_set = *Math::GSL::QRNGc::gsl_qrng_type_max_dimension_set;
*swig_state_size_get = *Math::GSL::QRNGc::gsl_qrng_type_state_size_get;
*swig_state_size_set = *Math::GSL::QRNGc::gsl_qrng_type_state_size_set;
*swig_init_state_get = *Math::GSL::QRNGc::gsl_qrng_type_init_state_get;
*swig_init_state_set = *Math::GSL::QRNGc::gsl_qrng_type_init_state_set;
*swig_get_get = *Math::GSL::QRNGc::gsl_qrng_type_get_get;
*swig_get_set = *Math::GSL::QRNGc::gsl_qrng_type_get_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::QRNGc::new_gsl_qrng_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::QRNGc::delete_gsl_qrng_type($self);
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


############# Class : Math::GSL::QRNG::gsl_qrng ##############

package Math::GSL::QRNG::gsl_qrng;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::QRNG );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::QRNGc::gsl_qrng_type_get;
*swig_type_set = *Math::GSL::QRNGc::gsl_qrng_type_set;
*swig_dimension_get = *Math::GSL::QRNGc::gsl_qrng_dimension_get;
*swig_dimension_set = *Math::GSL::QRNGc::gsl_qrng_dimension_set;
*swig_state_size_get = *Math::GSL::QRNGc::gsl_qrng_state_size_get;
*swig_state_size_set = *Math::GSL::QRNGc::gsl_qrng_state_size_set;
*swig_state_get = *Math::GSL::QRNGc::gsl_qrng_state_get;
*swig_state_set = *Math::GSL::QRNGc::gsl_qrng_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::QRNGc::new_gsl_qrng(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::QRNGc::delete_gsl_qrng($self);
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

package Math::GSL::QRNG;

*GSL_MAJOR_VERSION = *Math::GSL::QRNGc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::QRNGc::GSL_MINOR_VERSION;

my %__gsl_qrng_niederreiter_2_hash;
tie %__gsl_qrng_niederreiter_2_hash,"Math::GSL::QRNG::gsl_qrng_type", $Math::GSL::QRNGc::gsl_qrng_niederreiter_2;
$gsl_qrng_niederreiter_2= \%__gsl_qrng_niederreiter_2_hash;
bless $gsl_qrng_niederreiter_2, Math::GSL::QRNG::gsl_qrng_type;

my %__gsl_qrng_sobol_hash;
tie %__gsl_qrng_sobol_hash,"Math::GSL::QRNG::gsl_qrng_type", $Math::GSL::QRNGc::gsl_qrng_sobol;
$gsl_qrng_sobol= \%__gsl_qrng_sobol_hash;
bless $gsl_qrng_sobol, Math::GSL::QRNG::gsl_qrng_type;

my %__gsl_qrng_halton_hash;
tie %__gsl_qrng_halton_hash,"Math::GSL::QRNG::gsl_qrng_type", $Math::GSL::QRNGc::gsl_qrng_halton;
$gsl_qrng_halton= \%__gsl_qrng_halton_hash;
bless $gsl_qrng_halton, Math::GSL::QRNG::gsl_qrng_type;

my %__gsl_qrng_reversehalton_hash;
tie %__gsl_qrng_reversehalton_hash,"Math::GSL::QRNG::gsl_qrng_type", $Math::GSL::QRNGc::gsl_qrng_reversehalton;
$gsl_qrng_reversehalton= \%__gsl_qrng_reversehalton_hash;
bless $gsl_qrng_reversehalton, Math::GSL::QRNG::gsl_qrng_type;


@EXPORT_OK = qw($gsl_qrng_niederreiter_2 $gsl_qrng_sobol $gsl_qrng_halton $gsl_qrng_reversehalton
                gsl_qrng_alloc gsl_qrng_memcpy gsl_qrng_clone
                gsl_qrng_free  gsl_qrng_init gsl_qrng_name 
                gsl_qrng_size gsl_qrng_state gsl_qrng_get
            );
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );


__END__

=head1 NAME

Math::GSL::QRNG - Quasi-random number generator

=head1 SYNOPSIS

use Math::GSL::QRNG qw/:all/;

=head1 DESCRIPTION

Here is a list of all the functions included in this module :

=over

=item C<gsl_qrng_alloc($T, $n)> - This function returns a pointer to a newly-created instance of a quasi-random sequence generator of type $T and dimension $d. The type $T must be one of the constants included in this module.

=item C<gsl_qrng_clone($q)> - This function returns a pointer to a newly created generator which is an exact copy of the generator $q.

=item C<gsl_qrng_memcpy($dest, $src)> - This function copies the quasi-random sequence generator $src into the pre-existing generator $dest, making $dest into an exact copy of $src. The two generators must be of the same type.

=item C<gsl_qrng_free($q)> - This function frees all the memory associated with the generator $q. 

=item C<gsl_qrng_init($q)> - This function reinitializes the generator $q to its starting point. Note that quasi-random sequences do not use a seed and always produce the same set of values. 

=item C<gsl_qrng_name($q)> - This function returns a pointer to the name of the generator $q. 

=item C<gsl_qrng_size($q)> - This function returns the size of the state of generator r from the generator $q. You can use this information to access the state directly.

=item C<gsl_qrng_state($q)> - This function returns a pointer to the state of generator r from the generator $q. You can use this information to access the state directly.

=item C<gsl_qrng_get>

=back

This module also contains the following constants : 

=over

=item C<$gsl_qrng_niederreiter_2>

=item C<$gsl_qrng_sobol> 

=item C<$gsl_qrng_halton> 

=item C<$gsl_qrng_reversehalton>

=back

For more informations on the functions, we refer you to the GSL offcial documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want


=head1 EXAMPLES

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
1;

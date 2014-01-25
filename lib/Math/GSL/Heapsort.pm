# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Heapsort;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Heapsortc;
bootstrap Math::GSL::Heapsort;
package Math::GSL::Heapsort;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Heapsort;

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

package Math::GSL::Heapsort;

*gsl_permutation_alloc = *Math::GSL::Heapsortc::gsl_permutation_alloc;
*gsl_permutation_calloc = *Math::GSL::Heapsortc::gsl_permutation_calloc;
*gsl_permutation_init = *Math::GSL::Heapsortc::gsl_permutation_init;
*gsl_permutation_free = *Math::GSL::Heapsortc::gsl_permutation_free;
*gsl_permutation_memcpy = *Math::GSL::Heapsortc::gsl_permutation_memcpy;
*gsl_permutation_fread = *Math::GSL::Heapsortc::gsl_permutation_fread;
*gsl_permutation_fwrite = *Math::GSL::Heapsortc::gsl_permutation_fwrite;
*gsl_permutation_fscanf = *Math::GSL::Heapsortc::gsl_permutation_fscanf;
*gsl_permutation_fprintf = *Math::GSL::Heapsortc::gsl_permutation_fprintf;
*gsl_permutation_size = *Math::GSL::Heapsortc::gsl_permutation_size;
*gsl_permutation_data = *Math::GSL::Heapsortc::gsl_permutation_data;
*gsl_permutation_swap = *Math::GSL::Heapsortc::gsl_permutation_swap;
*gsl_permutation_valid = *Math::GSL::Heapsortc::gsl_permutation_valid;
*gsl_permutation_reverse = *Math::GSL::Heapsortc::gsl_permutation_reverse;
*gsl_permutation_inverse = *Math::GSL::Heapsortc::gsl_permutation_inverse;
*gsl_permutation_next = *Math::GSL::Heapsortc::gsl_permutation_next;
*gsl_permutation_prev = *Math::GSL::Heapsortc::gsl_permutation_prev;
*gsl_permutation_mul = *Math::GSL::Heapsortc::gsl_permutation_mul;
*gsl_permutation_linear_to_canonical = *Math::GSL::Heapsortc::gsl_permutation_linear_to_canonical;
*gsl_permutation_canonical_to_linear = *Math::GSL::Heapsortc::gsl_permutation_canonical_to_linear;
*gsl_permutation_inversions = *Math::GSL::Heapsortc::gsl_permutation_inversions;
*gsl_permutation_linear_cycles = *Math::GSL::Heapsortc::gsl_permutation_linear_cycles;
*gsl_permutation_canonical_cycles = *Math::GSL::Heapsortc::gsl_permutation_canonical_cycles;
*gsl_permutation_get = *Math::GSL::Heapsortc::gsl_permutation_get;
*gsl_heapsort = *Math::GSL::Heapsortc::gsl_heapsort;
*gsl_heapsort_index = *Math::GSL::Heapsortc::gsl_heapsort_index;

############# Class : Math::GSL::Heapsort::gsl_permutation_struct ##############

package Math::GSL::Heapsort::gsl_permutation_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Heapsort );
%OWNER = ();
%ITERATORS = ();
*swig_size_get = *Math::GSL::Heapsortc::gsl_permutation_struct_size_get;
*swig_size_set = *Math::GSL::Heapsortc::gsl_permutation_struct_size_set;
*swig_data_get = *Math::GSL::Heapsortc::gsl_permutation_struct_data_get;
*swig_data_set = *Math::GSL::Heapsortc::gsl_permutation_struct_data_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Heapsortc::new_gsl_permutation_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Heapsortc::delete_gsl_permutation_struct($self);
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

package Math::GSL::Heapsort;

*GSL_MAJOR_VERSION = *Math::GSL::Heapsortc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Heapsortc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Heapsortc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Heapsortc::GSL_NEGZERO;

@EXPORT_OK = qw/
               gsl_heapsort 
               gsl_heapsort_index 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Heapsort - Functions for sorting data, both directly and indirectly (using an index)

=head1 SYNOPSIS

This module is not yet implemented. Patches Welcome!

    use Math::GSL::Heapsort qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over

=item * gsl_heapsort 

=item * gsl_heapsort_index 

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

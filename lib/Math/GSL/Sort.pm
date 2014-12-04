# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Sort;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Sortc;
bootstrap Math::GSL::Sort;
package Math::GSL::Sort;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Sort;

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

package Math::GSL::Sort;

*gsl_sort = *Math::GSL::Sortc::gsl_sort;
*gsl_sort_index = *Math::GSL::Sortc::gsl_sort_index;
*gsl_sort_smallest = *Math::GSL::Sortc::gsl_sort_smallest;
*gsl_sort_smallest_index = *Math::GSL::Sortc::gsl_sort_smallest_index;
*gsl_sort_largest = *Math::GSL::Sortc::gsl_sort_largest;
*gsl_sort_largest_index = *Math::GSL::Sortc::gsl_sort_largest_index;
*gsl_sort_int = *Math::GSL::Sortc::gsl_sort_int;
*gsl_sort_int_index = *Math::GSL::Sortc::gsl_sort_int_index;
*gsl_sort_int_smallest = *Math::GSL::Sortc::gsl_sort_int_smallest;
*gsl_sort_int_smallest_index = *Math::GSL::Sortc::gsl_sort_int_smallest_index;
*gsl_sort_int_largest = *Math::GSL::Sortc::gsl_sort_int_largest;
*gsl_sort_int_largest_index = *Math::GSL::Sortc::gsl_sort_int_largest_index;
*gsl_sort_vector = *Math::GSL::Sortc::gsl_sort_vector;
*gsl_sort_vector2 = *Math::GSL::Sortc::gsl_sort_vector2;
*gsl_sort_vector_index = *Math::GSL::Sortc::gsl_sort_vector_index;
*gsl_sort_vector_smallest = *Math::GSL::Sortc::gsl_sort_vector_smallest;
*gsl_sort_vector_largest = *Math::GSL::Sortc::gsl_sort_vector_largest;
*gsl_sort_vector_smallest_index = *Math::GSL::Sortc::gsl_sort_vector_smallest_index;
*gsl_sort_vector_largest_index = *Math::GSL::Sortc::gsl_sort_vector_largest_index;
*gsl_sort_vector_int = *Math::GSL::Sortc::gsl_sort_vector_int;
*gsl_sort_vector_int_index = *Math::GSL::Sortc::gsl_sort_vector_int_index;
*gsl_sort_vector_int_smallest = *Math::GSL::Sortc::gsl_sort_vector_int_smallest;
*gsl_sort_vector_int_largest = *Math::GSL::Sortc::gsl_sort_vector_int_largest;
*gsl_sort_vector_int_smallest_index = *Math::GSL::Sortc::gsl_sort_vector_int_smallest_index;
*gsl_sort_vector_int_largest_index = *Math::GSL::Sortc::gsl_sort_vector_int_largest_index;
*gsl_permutation_alloc = *Math::GSL::Sortc::gsl_permutation_alloc;
*gsl_permutation_calloc = *Math::GSL::Sortc::gsl_permutation_calloc;
*gsl_permutation_init = *Math::GSL::Sortc::gsl_permutation_init;
*gsl_permutation_free = *Math::GSL::Sortc::gsl_permutation_free;
*gsl_permutation_memcpy = *Math::GSL::Sortc::gsl_permutation_memcpy;
*gsl_permutation_fread = *Math::GSL::Sortc::gsl_permutation_fread;
*gsl_permutation_fwrite = *Math::GSL::Sortc::gsl_permutation_fwrite;
*gsl_permutation_fscanf = *Math::GSL::Sortc::gsl_permutation_fscanf;
*gsl_permutation_fprintf = *Math::GSL::Sortc::gsl_permutation_fprintf;
*gsl_permutation_size = *Math::GSL::Sortc::gsl_permutation_size;
*gsl_permutation_data = *Math::GSL::Sortc::gsl_permutation_data;
*gsl_permutation_swap = *Math::GSL::Sortc::gsl_permutation_swap;
*gsl_permutation_valid = *Math::GSL::Sortc::gsl_permutation_valid;
*gsl_permutation_reverse = *Math::GSL::Sortc::gsl_permutation_reverse;
*gsl_permutation_inverse = *Math::GSL::Sortc::gsl_permutation_inverse;
*gsl_permutation_next = *Math::GSL::Sortc::gsl_permutation_next;
*gsl_permutation_prev = *Math::GSL::Sortc::gsl_permutation_prev;
*gsl_permutation_mul = *Math::GSL::Sortc::gsl_permutation_mul;
*gsl_permutation_linear_to_canonical = *Math::GSL::Sortc::gsl_permutation_linear_to_canonical;
*gsl_permutation_canonical_to_linear = *Math::GSL::Sortc::gsl_permutation_canonical_to_linear;
*gsl_permutation_inversions = *Math::GSL::Sortc::gsl_permutation_inversions;
*gsl_permutation_linear_cycles = *Math::GSL::Sortc::gsl_permutation_linear_cycles;
*gsl_permutation_canonical_cycles = *Math::GSL::Sortc::gsl_permutation_canonical_cycles;
*gsl_permutation_get = *Math::GSL::Sortc::gsl_permutation_get;

############# Class : Math::GSL::Sort::gsl_permutation_struct ##############

package Math::GSL::Sort::gsl_permutation_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Sort );
%OWNER = ();
%ITERATORS = ();
*swig_size_get = *Math::GSL::Sortc::gsl_permutation_struct_size_get;
*swig_size_set = *Math::GSL::Sortc::gsl_permutation_struct_size_set;
*swig_data_get = *Math::GSL::Sortc::gsl_permutation_struct_data_get;
*swig_data_set = *Math::GSL::Sortc::gsl_permutation_struct_data_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Sortc::new_gsl_permutation_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Sortc::delete_gsl_permutation_struct($self);
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

package Math::GSL::Sort;

*GSL_MAJOR_VERSION = *Math::GSL::Sortc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Sortc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Sortc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Sortc::GSL_NEGZERO;

@EXPORT_plain = qw/
                gsl_sort gsl_sort_index 
                gsl_sort_smallest gsl_sort_smallest_index
                gsl_sort_largest gsl_sort_largest_index
                /;
@EXPORT_vector= qw/
                gsl_sort_vector gsl_sort_vector_index 
                gsl_sort_vector_smallest gsl_sort_vector_smallest_index
                gsl_sort_vector_largest gsl_sort_vector_largest_index
                /;
@EXPORT_OK    = ( @EXPORT_plain, @EXPORT_vector );
%EXPORT_TAGS  = (
                 all    => [ @EXPORT_OK     ], 
                 plain  => [ @EXPORT_plain  ], 
                 vector => [ @EXPORT_vector ], 
                );
__END__

=encoding utf8

=head1 NAME

Math::GSL::Sort - Functions for sorting data

=head1 SYNOPSIS

    use Math::GSL::Sort qw/:all/;
    my $x       = [ 2**15, 1.67, 20e5, -17, 6900, 1/3 , 42e-10 ];
    my $sorted  = gsl_sort($x, 1, $#$x+1 );
    my $numbers = [ map { rand(100) } (1..100) ];
    my ($status, $smallest10) = gsl_sort_smallest($array, 10, $x, 1, $#$x+1);


=head1 DESCRIPTION

=over

=item * gsl_sort_vector($v) 

This function sorts the elements of the vector $v into ascending numerical
order.

=item * gsl_sort_vector_index($p, $v) 

This function indirectly sorts the elements of the vector $v into ascending
order, storing the resulting permutation in $p. The elements of $p give the
index of the vector element which would have been stored in that position if
the vector had been sorted in place. The first element of $p gives the index
of the least element in $v, and the last element of $p gives the index of the
greatest element in $v. The vector $v is not changed. 

=item * gsl_sort_vector_smallest($array, $k, $vector) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
$k smallest elements of the vector $v. $k must be less than or equal to the
length of the vector $v.

=item * gsl_sort_vector_smallest_index($p, $k, $v) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
indices of the $k smallest elements of the vector $v. $p must be a prealocated
array reference. This should be removed in further versions. $k must be less
than or equal to the length of the vector $v. 

=item * gsl_sort_vector_largest($array, $k, $vector) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
$k largest elements of the vector $v. $k must be less than or equal to the
length of the vector $v.

=item * gsl_sort_vector_largest_index($p, $k, $v) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
indices of the $k largest elements of the vector $v. $p must be a prealocated
array reference. This should be removed in further versions. $k must be less
than or equal to the length of the vector $v. 

=item * gsl_sort($data, $stride, $n) 

This function returns an array reference to the sorted $n elements of the
array $data with stride $stride into ascending numerical order.

=item * gsl_sort_index($p, $data, $stride, $n) 

This function indirectly sorts the $n elements of the array $data with stride
$stride into ascending order, outputting the permutation in the foram of an
array. $p must be a prealocated array reference. This should be removed in
further versions. The array $data is not changed.

=item * gsl_sort_smallest($array, $k, $data, $stride, $n) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
$k smallest elements of the array $data, of size $n and stride $stride, in
ascending numerical. The size $k of the subset must be less than or equal to
$n. The data $src is not modified by this operation. $array must be a
prealocated array reference. This should be removed in further versions.

=item * gsl_sort_smallest_index($p, $k, $src, $stride, $n) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
indices of the $k smallest elements of the array $src, of size $n and stride
$stride. The indices are chosen so that the corresponding data is in ascending
numerical order. $k must be less than or equal to $n. The data $src is not
modified by this operation. $p must be a prealocated array reference. This
should be removed in further versions. 

=item * gsl_sort_largest($array, $k, $data, $stride, $n) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
$k largest elements of the array $data, of size $n and stride $stride, in
ascending numerical. The size $k of the subset must be less than or equal to
$n. The data $src is not modified by this operation. $array must be a
prealocated array reference. This should be removed in further versions.

=item * gsl_sort_largest_index($p, $k, $src, $stride, $n) 

This function outputs 0 if the operation succeeded, 1 otherwise and then the
indices of the $k largest elements of the array $src, of size $n and stride
$stride. The indices are chosen so that the corresponding data is in ascending
numerical order. $k must be less than or equal to $n. The data $src is not
modified by this operation. $p must be a prealocated array reference. This
should be removed in further versions. 

=back

 Here is a complete list of all tags for this module :

=over

=item all

=item plain

=item vector

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

=head1 PERFORMANCE

In the source code of Math::GSL, the file "examples/benchmark/sort" compares
the performance of gsl_sort() to Perl's builtin sort() function. It's first
argument is the number of iterations and the second is the size of the array
of numbers to sort. For example, to see a benchmark of 1000 iterations for 
arrays of size 50000 you would type

    ./examples/benchmark/sort 1000 50000

Initial benchmarks indicate just slightly above a 2x performance increase
over sort() for arrays of between 5000 and 50000 elements. This may mostly
be due to the fact that gsl_sort() takes and returns a reference while sort()
takes and returns a plain list.

=head1 AUTHORS

Jonathan "Duke" Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2011 Jonathan "Duke" Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::NTuple;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::NTuplec;
bootstrap Math::GSL::NTuple;
package Math::GSL::NTuple;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::NTuple;

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

package Math::GSL::NTuple;

*gsl_ntuple_open = *Math::GSL::NTuplec::gsl_ntuple_open;
*gsl_ntuple_create = *Math::GSL::NTuplec::gsl_ntuple_create;
*gsl_ntuple_write = *Math::GSL::NTuplec::gsl_ntuple_write;
*gsl_ntuple_read = *Math::GSL::NTuplec::gsl_ntuple_read;
*gsl_ntuple_bookdata = *Math::GSL::NTuplec::gsl_ntuple_bookdata;
*gsl_ntuple_project = *Math::GSL::NTuplec::gsl_ntuple_project;
*gsl_ntuple_close = *Math::GSL::NTuplec::gsl_ntuple_close;
*gsl_histogram_alloc = *Math::GSL::NTuplec::gsl_histogram_alloc;
*gsl_histogram_calloc = *Math::GSL::NTuplec::gsl_histogram_calloc;
*gsl_histogram_calloc_uniform = *Math::GSL::NTuplec::gsl_histogram_calloc_uniform;
*gsl_histogram_free = *Math::GSL::NTuplec::gsl_histogram_free;
*gsl_histogram_increment = *Math::GSL::NTuplec::gsl_histogram_increment;
*gsl_histogram_accumulate = *Math::GSL::NTuplec::gsl_histogram_accumulate;
*gsl_histogram_find = *Math::GSL::NTuplec::gsl_histogram_find;
*gsl_histogram_get = *Math::GSL::NTuplec::gsl_histogram_get;
*gsl_histogram_get_range = *Math::GSL::NTuplec::gsl_histogram_get_range;
*gsl_histogram_max = *Math::GSL::NTuplec::gsl_histogram_max;
*gsl_histogram_min = *Math::GSL::NTuplec::gsl_histogram_min;
*gsl_histogram_bins = *Math::GSL::NTuplec::gsl_histogram_bins;
*gsl_histogram_reset = *Math::GSL::NTuplec::gsl_histogram_reset;
*gsl_histogram_calloc_range = *Math::GSL::NTuplec::gsl_histogram_calloc_range;
*gsl_histogram_set_ranges = *Math::GSL::NTuplec::gsl_histogram_set_ranges;
*gsl_histogram_set_ranges_uniform = *Math::GSL::NTuplec::gsl_histogram_set_ranges_uniform;
*gsl_histogram_memcpy = *Math::GSL::NTuplec::gsl_histogram_memcpy;
*gsl_histogram_clone = *Math::GSL::NTuplec::gsl_histogram_clone;
*gsl_histogram_max_val = *Math::GSL::NTuplec::gsl_histogram_max_val;
*gsl_histogram_max_bin = *Math::GSL::NTuplec::gsl_histogram_max_bin;
*gsl_histogram_min_val = *Math::GSL::NTuplec::gsl_histogram_min_val;
*gsl_histogram_min_bin = *Math::GSL::NTuplec::gsl_histogram_min_bin;
*gsl_histogram_equal_bins_p = *Math::GSL::NTuplec::gsl_histogram_equal_bins_p;
*gsl_histogram_add = *Math::GSL::NTuplec::gsl_histogram_add;
*gsl_histogram_sub = *Math::GSL::NTuplec::gsl_histogram_sub;
*gsl_histogram_mul = *Math::GSL::NTuplec::gsl_histogram_mul;
*gsl_histogram_div = *Math::GSL::NTuplec::gsl_histogram_div;
*gsl_histogram_scale = *Math::GSL::NTuplec::gsl_histogram_scale;
*gsl_histogram_shift = *Math::GSL::NTuplec::gsl_histogram_shift;
*gsl_histogram_sigma = *Math::GSL::NTuplec::gsl_histogram_sigma;
*gsl_histogram_mean = *Math::GSL::NTuplec::gsl_histogram_mean;
*gsl_histogram_sum = *Math::GSL::NTuplec::gsl_histogram_sum;
*gsl_histogram_fwrite = *Math::GSL::NTuplec::gsl_histogram_fwrite;
*gsl_histogram_fread = *Math::GSL::NTuplec::gsl_histogram_fread;
*gsl_histogram_fprintf = *Math::GSL::NTuplec::gsl_histogram_fprintf;
*gsl_histogram_fscanf = *Math::GSL::NTuplec::gsl_histogram_fscanf;
*gsl_histogram_pdf_alloc = *Math::GSL::NTuplec::gsl_histogram_pdf_alloc;
*gsl_histogram_pdf_init = *Math::GSL::NTuplec::gsl_histogram_pdf_init;
*gsl_histogram_pdf_free = *Math::GSL::NTuplec::gsl_histogram_pdf_free;
*gsl_histogram_pdf_sample = *Math::GSL::NTuplec::gsl_histogram_pdf_sample;

############# Class : Math::GSL::NTuple::gsl_ntuple ##############

package Math::GSL::NTuple::gsl_ntuple;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::NTuple );
%OWNER = ();
%ITERATORS = ();
*swig_file_get = *Math::GSL::NTuplec::gsl_ntuple_file_get;
*swig_file_set = *Math::GSL::NTuplec::gsl_ntuple_file_set;
*swig_ntuple_data_get = *Math::GSL::NTuplec::gsl_ntuple_ntuple_data_get;
*swig_ntuple_data_set = *Math::GSL::NTuplec::gsl_ntuple_ntuple_data_set;
*swig_size_get = *Math::GSL::NTuplec::gsl_ntuple_size_get;
*swig_size_set = *Math::GSL::NTuplec::gsl_ntuple_size_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::NTuplec::new_gsl_ntuple(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::NTuplec::delete_gsl_ntuple($self);
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


############# Class : Math::GSL::NTuple::gsl_ntuple_select_fn ##############

package Math::GSL::NTuple::gsl_ntuple_select_fn;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::NTuple );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::NTuplec::gsl_ntuple_select_fn_function_get;
*swig_function_set = *Math::GSL::NTuplec::gsl_ntuple_select_fn_function_set;
*swig_params_get = *Math::GSL::NTuplec::gsl_ntuple_select_fn_params_get;
*swig_params_set = *Math::GSL::NTuplec::gsl_ntuple_select_fn_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::NTuplec::new_gsl_ntuple_select_fn(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::NTuplec::delete_gsl_ntuple_select_fn($self);
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


############# Class : Math::GSL::NTuple::gsl_ntuple_value_fn ##############

package Math::GSL::NTuple::gsl_ntuple_value_fn;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::NTuple );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::NTuplec::gsl_ntuple_value_fn_function_get;
*swig_function_set = *Math::GSL::NTuplec::gsl_ntuple_value_fn_function_set;
*swig_params_get = *Math::GSL::NTuplec::gsl_ntuple_value_fn_params_get;
*swig_params_set = *Math::GSL::NTuplec::gsl_ntuple_value_fn_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::NTuplec::new_gsl_ntuple_value_fn(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::NTuplec::delete_gsl_ntuple_value_fn($self);
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


############# Class : Math::GSL::NTuple::gsl_histogram ##############

package Math::GSL::NTuple::gsl_histogram;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::NTuple );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::NTuplec::gsl_histogram_n_get;
*swig_n_set = *Math::GSL::NTuplec::gsl_histogram_n_set;
*swig_range_get = *Math::GSL::NTuplec::gsl_histogram_range_get;
*swig_range_set = *Math::GSL::NTuplec::gsl_histogram_range_set;
*swig_bin_get = *Math::GSL::NTuplec::gsl_histogram_bin_get;
*swig_bin_set = *Math::GSL::NTuplec::gsl_histogram_bin_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::NTuplec::new_gsl_histogram(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::NTuplec::delete_gsl_histogram($self);
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


############# Class : Math::GSL::NTuple::gsl_histogram_pdf ##############

package Math::GSL::NTuple::gsl_histogram_pdf;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::NTuple );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::NTuplec::gsl_histogram_pdf_n_get;
*swig_n_set = *Math::GSL::NTuplec::gsl_histogram_pdf_n_set;
*swig_range_get = *Math::GSL::NTuplec::gsl_histogram_pdf_range_get;
*swig_range_set = *Math::GSL::NTuplec::gsl_histogram_pdf_range_set;
*swig_sum_get = *Math::GSL::NTuplec::gsl_histogram_pdf_sum_get;
*swig_sum_set = *Math::GSL::NTuplec::gsl_histogram_pdf_sum_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::NTuplec::new_gsl_histogram_pdf(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::NTuplec::delete_gsl_histogram_pdf($self);
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

package Math::GSL::NTuple;

*GSL_MAJOR_VERSION = *Math::GSL::NTuplec::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::NTuplec::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::NTuplec::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::NTuplec::GSL_NEGZERO;

use Data::Dumper;
use Carp qw/croak/;
use Math::GSL::Errno qw/:all/;

@EXPORT_OK = qw/
               gsl_ntuple_open
               gsl_ntuple_create
               gsl_ntuple_write
               gsl_ntuple_read
               gsl_ntuple_bookdata
               gsl_ntuple_project
               gsl_ntuple_close
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

=encoding utf8

=head1 NAME

Math::GSL::NTuple - Functions for creating and manipulating ntuples, sets of values

=head1 SYNOPSIS

This module is partially implemented. Patches Welcome!

    use Math::GSL::NTuple qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over

=cut

sub new
{
    my ($class,$values) = @_;
    my $this = {};
    my $ntuple = Math::GSL::NTuple::gsl_ntuple->new;
    $this->{_ntuple} = $ntuple;

    bless $this, $class;
}

sub raw
{
    return (shift)->{_ntuple};
}

=item * <gsl_ntuple_open($filename, $ntuple_data, $size)>

This function opens an existing ntuple file $filename for reading and returns a
pointer to a corresponding ntuple struct. The ntuples in the file must have size
$size. A pointer to memory for the current ntuple row $ntuple_data, which is an
array reference, must be supplied -- this is used to copy ntuples in and out of
the file.

=item * <gsl_ntuple_create>

This function creates a new write-only ntuple file $filename for ntuples of size
$size and returns a pointer to the newly created ntuple struct. Any existing
file with the same name is truncated to zero length and overwritten. A pointer
to memory for the current ntuple row $ntuple_data, which is an array reference,
must be supplied -- this is used to copy ntuples in and out of the file.


=item * <gsl_ntuple_write($ntuple)>

This function writes the current $ntuple $ntuple->{ntuple_data} of size
$ntuple->{size} to the corresponding file.

=item * <gsl_ntuple_bookdata($ntuple)>

This function is a synonym for gsl_ntuple_write.

=item * <gsl_ntuple_read($ntuple)>

This function reads the current row of the ntuple file for ntuple and stores the
values in $ntuple->{data}.

=item * <gsl_ntuple_project()>

=item * <gsl_ntuple_close($ntuple)> 

This function closes the ntuple file ntuple and frees its associated allocated
memory.

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

=head1 AUTHORS

Jonathan "Duke" Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2014 Jonathan "Duke" Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

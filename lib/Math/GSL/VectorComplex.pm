# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::VectorComplex;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::VectorComplexc;
bootstrap Math::GSL::VectorComplex;
package Math::GSL::VectorComplex;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::VectorComplex;

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

package Math::GSL::VectorComplex;

*gsl_vector_alloc = *Math::GSL::VectorComplexc::gsl_vector_alloc;
*gsl_vector_calloc = *Math::GSL::VectorComplexc::gsl_vector_calloc;
*gsl_vector_alloc_from_block = *Math::GSL::VectorComplexc::gsl_vector_alloc_from_block;
*gsl_vector_alloc_from_vector = *Math::GSL::VectorComplexc::gsl_vector_alloc_from_vector;
*gsl_vector_free = *Math::GSL::VectorComplexc::gsl_vector_free;
*gsl_vector_view_array = *Math::GSL::VectorComplexc::gsl_vector_view_array;
*gsl_vector_view_array_with_stride = *Math::GSL::VectorComplexc::gsl_vector_view_array_with_stride;
*gsl_vector_const_view_array = *Math::GSL::VectorComplexc::gsl_vector_const_view_array;
*gsl_vector_const_view_array_with_stride = *Math::GSL::VectorComplexc::gsl_vector_const_view_array_with_stride;
*gsl_vector_subvector = *Math::GSL::VectorComplexc::gsl_vector_subvector;
*gsl_vector_subvector_with_stride = *Math::GSL::VectorComplexc::gsl_vector_subvector_with_stride;
*gsl_vector_const_subvector = *Math::GSL::VectorComplexc::gsl_vector_const_subvector;
*gsl_vector_const_subvector_with_stride = *Math::GSL::VectorComplexc::gsl_vector_const_subvector_with_stride;
*gsl_vector_set_zero = *Math::GSL::VectorComplexc::gsl_vector_set_zero;
*gsl_vector_set_all = *Math::GSL::VectorComplexc::gsl_vector_set_all;
*gsl_vector_set_basis = *Math::GSL::VectorComplexc::gsl_vector_set_basis;
*gsl_vector_fread = *Math::GSL::VectorComplexc::gsl_vector_fread;
*gsl_vector_fwrite = *Math::GSL::VectorComplexc::gsl_vector_fwrite;
*gsl_vector_fscanf = *Math::GSL::VectorComplexc::gsl_vector_fscanf;
*gsl_vector_fprintf = *Math::GSL::VectorComplexc::gsl_vector_fprintf;
*gsl_vector_memcpy = *Math::GSL::VectorComplexc::gsl_vector_memcpy;
*gsl_vector_reverse = *Math::GSL::VectorComplexc::gsl_vector_reverse;
*gsl_vector_swap = *Math::GSL::VectorComplexc::gsl_vector_swap;
*gsl_vector_swap_elements = *Math::GSL::VectorComplexc::gsl_vector_swap_elements;
*gsl_vector_max = *Math::GSL::VectorComplexc::gsl_vector_max;
*gsl_vector_min = *Math::GSL::VectorComplexc::gsl_vector_min;
*gsl_vector_minmax = *Math::GSL::VectorComplexc::gsl_vector_minmax;
*gsl_vector_max_index = *Math::GSL::VectorComplexc::gsl_vector_max_index;
*gsl_vector_min_index = *Math::GSL::VectorComplexc::gsl_vector_min_index;
*gsl_vector_minmax_index = *Math::GSL::VectorComplexc::gsl_vector_minmax_index;
*gsl_vector_add = *Math::GSL::VectorComplexc::gsl_vector_add;
*gsl_vector_sub = *Math::GSL::VectorComplexc::gsl_vector_sub;
*gsl_vector_mul = *Math::GSL::VectorComplexc::gsl_vector_mul;
*gsl_vector_div = *Math::GSL::VectorComplexc::gsl_vector_div;
*gsl_vector_scale = *Math::GSL::VectorComplexc::gsl_vector_scale;
*gsl_vector_add_constant = *Math::GSL::VectorComplexc::gsl_vector_add_constant;
*gsl_vector_isnull = *Math::GSL::VectorComplexc::gsl_vector_isnull;
*gsl_vector_ispos = *Math::GSL::VectorComplexc::gsl_vector_ispos;
*gsl_vector_isneg = *Math::GSL::VectorComplexc::gsl_vector_isneg;
*gsl_vector_isnonneg = *Math::GSL::VectorComplexc::gsl_vector_isnonneg;
*gsl_vector_get = *Math::GSL::VectorComplexc::gsl_vector_get;
*gsl_vector_set = *Math::GSL::VectorComplexc::gsl_vector_set;
*gsl_vector_ptr = *Math::GSL::VectorComplexc::gsl_vector_ptr;
*gsl_vector_const_ptr = *Math::GSL::VectorComplexc::gsl_vector_const_ptr;
*gsl_vector_complex_alloc = *Math::GSL::VectorComplexc::gsl_vector_complex_alloc;
*gsl_vector_complex_calloc = *Math::GSL::VectorComplexc::gsl_vector_complex_calloc;
*gsl_vector_complex_alloc_from_block = *Math::GSL::VectorComplexc::gsl_vector_complex_alloc_from_block;
*gsl_vector_complex_alloc_from_vector = *Math::GSL::VectorComplexc::gsl_vector_complex_alloc_from_vector;
*gsl_vector_complex_free = *Math::GSL::VectorComplexc::gsl_vector_complex_free;
*gsl_vector_complex_view_array = *Math::GSL::VectorComplexc::gsl_vector_complex_view_array;
*gsl_vector_complex_view_array_with_stride = *Math::GSL::VectorComplexc::gsl_vector_complex_view_array_with_stride;
*gsl_vector_complex_const_view_array = *Math::GSL::VectorComplexc::gsl_vector_complex_const_view_array;
*gsl_vector_complex_const_view_array_with_stride = *Math::GSL::VectorComplexc::gsl_vector_complex_const_view_array_with_stride;
*gsl_vector_complex_subvector = *Math::GSL::VectorComplexc::gsl_vector_complex_subvector;
*gsl_vector_complex_subvector_with_stride = *Math::GSL::VectorComplexc::gsl_vector_complex_subvector_with_stride;
*gsl_vector_complex_const_subvector = *Math::GSL::VectorComplexc::gsl_vector_complex_const_subvector;
*gsl_vector_complex_const_subvector_with_stride = *Math::GSL::VectorComplexc::gsl_vector_complex_const_subvector_with_stride;
*gsl_vector_complex_real = *Math::GSL::VectorComplexc::gsl_vector_complex_real;
*gsl_vector_complex_imag = *Math::GSL::VectorComplexc::gsl_vector_complex_imag;
*gsl_vector_complex_const_real = *Math::GSL::VectorComplexc::gsl_vector_complex_const_real;
*gsl_vector_complex_const_imag = *Math::GSL::VectorComplexc::gsl_vector_complex_const_imag;
*gsl_vector_complex_set_zero = *Math::GSL::VectorComplexc::gsl_vector_complex_set_zero;
*gsl_vector_complex_set_all = *Math::GSL::VectorComplexc::gsl_vector_complex_set_all;
*gsl_vector_complex_set_basis = *Math::GSL::VectorComplexc::gsl_vector_complex_set_basis;
*gsl_vector_complex_fread = *Math::GSL::VectorComplexc::gsl_vector_complex_fread;
*gsl_vector_complex_fwrite = *Math::GSL::VectorComplexc::gsl_vector_complex_fwrite;
*gsl_vector_complex_fscanf = *Math::GSL::VectorComplexc::gsl_vector_complex_fscanf;
*gsl_vector_complex_fprintf = *Math::GSL::VectorComplexc::gsl_vector_complex_fprintf;
*gsl_vector_complex_memcpy = *Math::GSL::VectorComplexc::gsl_vector_complex_memcpy;
*gsl_vector_complex_reverse = *Math::GSL::VectorComplexc::gsl_vector_complex_reverse;
*gsl_vector_complex_swap = *Math::GSL::VectorComplexc::gsl_vector_complex_swap;
*gsl_vector_complex_swap_elements = *Math::GSL::VectorComplexc::gsl_vector_complex_swap_elements;
*gsl_vector_complex_isnull = *Math::GSL::VectorComplexc::gsl_vector_complex_isnull;
*gsl_vector_complex_ispos = *Math::GSL::VectorComplexc::gsl_vector_complex_ispos;
*gsl_vector_complex_isneg = *Math::GSL::VectorComplexc::gsl_vector_complex_isneg;
*gsl_vector_complex_isnonneg = *Math::GSL::VectorComplexc::gsl_vector_complex_isnonneg;
*gsl_vector_complex_add = *Math::GSL::VectorComplexc::gsl_vector_complex_add;
*gsl_vector_complex_sub = *Math::GSL::VectorComplexc::gsl_vector_complex_sub;
*gsl_vector_complex_mul = *Math::GSL::VectorComplexc::gsl_vector_complex_mul;
*gsl_vector_complex_div = *Math::GSL::VectorComplexc::gsl_vector_complex_div;
*gsl_vector_complex_scale = *Math::GSL::VectorComplexc::gsl_vector_complex_scale;
*gsl_vector_complex_add_constant = *Math::GSL::VectorComplexc::gsl_vector_complex_add_constant;
*gsl_vector_complex_get = *Math::GSL::VectorComplexc::gsl_vector_complex_get;
*gsl_vector_complex_set = *Math::GSL::VectorComplexc::gsl_vector_complex_set;
*gsl_vector_complex_ptr = *Math::GSL::VectorComplexc::gsl_vector_complex_ptr;
*gsl_vector_complex_const_ptr = *Math::GSL::VectorComplexc::gsl_vector_complex_const_ptr;

############# Class : Math::GSL::VectorComplex::gsl_complex_long_double ##############

package Math::GSL::VectorComplex::gsl_complex_long_double;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::VectorComplexc::gsl_complex_long_double_dat_get;
*swig_dat_set = *Math::GSL::VectorComplexc::gsl_complex_long_double_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new_gsl_complex_long_double(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete_gsl_complex_long_double($self);
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


############# Class : Math::GSL::VectorComplex::gsl_complex ##############

package Math::GSL::VectorComplex::gsl_complex;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::VectorComplexc::gsl_complex_dat_get;
*swig_dat_set = *Math::GSL::VectorComplexc::gsl_complex_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new_gsl_complex(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete_gsl_complex($self);
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


############# Class : Math::GSL::VectorComplex::gsl_complex_float ##############

package Math::GSL::VectorComplex::gsl_complex_float;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::VectorComplexc::gsl_complex_float_dat_get;
*swig_dat_set = *Math::GSL::VectorComplexc::gsl_complex_float_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new_gsl_complex_float(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete_gsl_complex_float($self);
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


############# Class : Math::GSL::VectorComplex::gsl_vector ##############

package Math::GSL::VectorComplex::gsl_vector;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_size_get = *Math::GSL::VectorComplexc::gsl_vector_size_get;
*swig_size_set = *Math::GSL::VectorComplexc::gsl_vector_size_set;
*swig_stride_get = *Math::GSL::VectorComplexc::gsl_vector_stride_get;
*swig_stride_set = *Math::GSL::VectorComplexc::gsl_vector_stride_set;
*swig_data_get = *Math::GSL::VectorComplexc::gsl_vector_data_get;
*swig_data_set = *Math::GSL::VectorComplexc::gsl_vector_data_set;
*swig_block_get = *Math::GSL::VectorComplexc::gsl_vector_block_get;
*swig_block_set = *Math::GSL::VectorComplexc::gsl_vector_block_set;
*swig_owner_get = *Math::GSL::VectorComplexc::gsl_vector_owner_get;
*swig_owner_set = *Math::GSL::VectorComplexc::gsl_vector_owner_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new_gsl_vector(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete_gsl_vector($self);
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


############# Class : Math::GSL::VectorComplex::_gsl_vector_view ##############

package Math::GSL::VectorComplex::_gsl_vector_view;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_vector_get = *Math::GSL::VectorComplexc::_gsl_vector_view_vector_get;
*swig_vector_set = *Math::GSL::VectorComplexc::_gsl_vector_view_vector_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new__gsl_vector_view(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete__gsl_vector_view($self);
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


############# Class : Math::GSL::VectorComplex::_gsl_vector_const_view ##############

package Math::GSL::VectorComplex::_gsl_vector_const_view;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_vector_get = *Math::GSL::VectorComplexc::_gsl_vector_const_view_vector_get;
*swig_vector_set = *Math::GSL::VectorComplexc::_gsl_vector_const_view_vector_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new__gsl_vector_const_view(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete__gsl_vector_const_view($self);
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


############# Class : Math::GSL::VectorComplex::gsl_vector_complex ##############

package Math::GSL::VectorComplex::gsl_vector_complex;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_size_get = *Math::GSL::VectorComplexc::gsl_vector_complex_size_get;
*swig_size_set = *Math::GSL::VectorComplexc::gsl_vector_complex_size_set;
*swig_stride_get = *Math::GSL::VectorComplexc::gsl_vector_complex_stride_get;
*swig_stride_set = *Math::GSL::VectorComplexc::gsl_vector_complex_stride_set;
*swig_data_get = *Math::GSL::VectorComplexc::gsl_vector_complex_data_get;
*swig_data_set = *Math::GSL::VectorComplexc::gsl_vector_complex_data_set;
*swig_block_get = *Math::GSL::VectorComplexc::gsl_vector_complex_block_get;
*swig_block_set = *Math::GSL::VectorComplexc::gsl_vector_complex_block_set;
*swig_owner_get = *Math::GSL::VectorComplexc::gsl_vector_complex_owner_get;
*swig_owner_set = *Math::GSL::VectorComplexc::gsl_vector_complex_owner_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new_gsl_vector_complex(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete_gsl_vector_complex($self);
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


############# Class : Math::GSL::VectorComplex::_gsl_vector_complex_view ##############

package Math::GSL::VectorComplex::_gsl_vector_complex_view;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_vector_get = *Math::GSL::VectorComplexc::_gsl_vector_complex_view_vector_get;
*swig_vector_set = *Math::GSL::VectorComplexc::_gsl_vector_complex_view_vector_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new__gsl_vector_complex_view(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete__gsl_vector_complex_view($self);
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


############# Class : Math::GSL::VectorComplex::_gsl_vector_complex_const_view ##############

package Math::GSL::VectorComplex::_gsl_vector_complex_const_view;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::VectorComplex );
%OWNER = ();
%ITERATORS = ();
*swig_vector_get = *Math::GSL::VectorComplexc::_gsl_vector_complex_const_view_vector_get;
*swig_vector_set = *Math::GSL::VectorComplexc::_gsl_vector_complex_const_view_vector_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::VectorComplexc::new__gsl_vector_complex_const_view(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::VectorComplexc::delete__gsl_vector_complex_const_view($self);
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

package Math::GSL::VectorComplex;

*GSL_MAJOR_VERSION = *Math::GSL::VectorComplexc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::VectorComplexc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::VectorComplexc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::VectorComplexc::GSL_NEGZERO;

use Scalar::Util 'blessed';
use Data::Dumper;
use Carp qw/croak/;
use Math::GSL::Errno qw/:all/;
use Math::GSL::BLAS qw/gsl_blas_ddot/;
use Math::GSL::Complex qw/:all/;
use Math::GSL::Test qw/is_status_ok/;
use Math::Complex;
use overload
    '*'      => \&_multiplication,
    '+'      => \&_addition,
    '-'      => \&_subtract,
    fallback => 1,
;

@EXPORT_all  = qw/
                 gsl_vector_complex_alloc gsl_vector_complex_calloc gsl_vector_complex_alloc_from_block gsl_vector_complex_alloc_from_vector
                 gsl_vector_complex_free gsl_vector_complex_view_array gsl_vector_complex_view_array_with_stride gsl_vector_complex_const_view_array
                 gsl_vector_complex_const_view_array_with_stride gsl_vector_complex_subvector gsl_vector_complex_subvector_with_stride
                 gsl_vector_complex_const_subvector gsl_vector_complex_const_subvector_with_stride gsl_vector_complex_real gsl_vector_complex_imag
                 gsl_vector_complex_const_real gsl_vector_complex_const_imag gsl_vector_complex_get gsl_vector_complex_set
                 gsl_vector_complex_ptr gsl_vector_complex_const_ptr gsl_vector_complex_set_zero gsl_vector_complex_set_all
                 gsl_vector_complex_set_basis gsl_vector_complex_fread gsl_vector_complex_fwrite gsl_vector_complex_fscanf
                 gsl_vector_complex_fprintf gsl_vector_complex_memcpy gsl_vector_complex_reverse gsl_vector_complex_swap
                 gsl_vector_complex_swap_elements gsl_vector_complex_isnull gsl_vector_complex_ispos gsl_vector_complex_isneg
/;
@EXPORT_OK = (@EXPORT_all);
%EXPORT_TAGS = ( all => \@EXPORT_all );

=head1 NAME

Math::GSL::VectorComplex - Complex Vectors

=head1 SYNOPSIS

    use Math::GSL::VectorComplex qw/:all/;
    my $vec1 = Math::GSL::VectorComplex->new([1 + 2*i, 7*i, 5, -3 ]);
    my $vec2 = $vec1 * 5;
    my $vec3 = Math::GSL::Vector>new(10);   # 10 element zero vector 
    my $vec4 = $vec1 + $vec2;

    # set the element at index 1 to -i
    # and the element at index 3 to i
    $vec3->set([ 1, -i ], [ 9, i ]);

    my @vec = $vec2->as_list;               # return elements as Perl list

    my $dot_product = $vec1 * $vec2;
    my $length      = $vec2->length;
    my $first       = $vec1->get(0);


=cut

=head1 Objected Oriented Interface to GSL Math::GSL::VectorComplex

=head2 new()

Creates a new Vector of the given size.

    my $vector = Math::GSL::VectorComplex->new(3);

You can also create and set directly the values of the vector like this :

   my $vector = Math::GSL::VectorComplex->new([2,4,1]);

=cut

sub new {
    my ($class, $values) = @_;
    my $length  = $#$values;
    my $this = {};
    my $vector;

    # we expect $values to have Math::Complex objects
    @$values = map { gsl_complex_rect(Re($_), Im($_)) } @$values;

    if ( ref $values eq 'ARRAY' ){
        die __PACKAGE__.'::new($x) - $x must be a nonempty array reference' if $length == -1;
        $vector  = gsl_vector_complex_alloc($length+1);
        map { gsl_vector_complex_set($vector, $_, $values->[$_] ) }  (0 .. $length);
        $this->{_length} = $length+1;
    } elsif ( (int($values) == $values) && ($values > 0)) {
        $vector  = gsl_vector_complex_alloc($values);
        gsl_vector_complex_set_zero($vector);
        $this->{_length} = $values;
    } else {
        die __PACKAGE__.'::new($x) - $x must be an int or array reference';
    }
    $this->{_vector} = $vector;
    bless $this, $class;
}


=head2 raw()

Get the underlying GSL vector object created by SWIG, useful for using gsl_vector_* functions which do not have an OO counterpart.

    my $vector    = Math::GSL::VectorComplex->new(3);
    my $gsl_vector = $vector->raw;
    my $stuff      = gsl_vector_get($gsl_vector, 1);

=cut

sub raw { 
    my $self = shift;
    return $self->{_vector};
}

=head2 min()

Returns the minimum value contained in the vector.

   my $vector = Math::GSL::VectorComplex->new([2,4,1]);
   my $minimum = $vector->min;

=cut 

sub min {
    my $self=shift;
    return gsl_vector_min($self->raw);
}

=head2 max()

Returns the minimum value contained in the vector.

   my $vector = Math::GSL::VectorComplex->new([2,4,1]);
   my $maximum = $vector->max;

=cut 

sub max {
    my $self=shift;
    return gsl_vector_max($self->raw);
}

=head2 length()

Returns the number of elements contained in the vector.

   my $vector = Math::GSL::VectorComplex->new([2,4,1]);
   my $length = $vector->length;

=cut 

sub length { my $self=shift; $self->{_length} }

=head2  as_list() 

Gets the content of a Math::GSL::Vector object as a Perl list.

    my $vector = Math::GSL::VectorComplex->new(3);
    ...
    my @values = $vector->as_list;
=cut

sub as_list {
    my $self=shift;
    # this is wrong
    return map { cplxe( gsl_complex_abs($_), gsl_complex_arg($_) ) } $self->get( [ 0 .. $self->length - 1  ] );
}

=head2  get()

Gets the value of an of a Math::GSL::Vector object.

    my $vector = Math::GSL::VectorComplex->new(3);
    ...
    my @values = $vector->get(2);

You can also enter an array of indices to receive their corresponding values:

    my $vector = Math::GSL::VectorComplex->new(3);
    ...
    my @values = $vector->get([0,2]);

=cut

sub get {
    my ($self, $indices) = @_;
    return  map {  gsl_vector_complex_get($self->raw, $_ ) } @$indices ;
}

=head2 reverse()

Returns the a vector with the elements in reversed order.

    use Math::Complex;
    my $v1 = Math::GSL::VectorComplex->new([ 1, 2, 3*i]);
    my $v2 = $v1->reverse;

=cut

sub reverse {
    my $self = shift;
    my $copy = $self->copy();
    unless ( is_status_ok( gsl_vector_complex_reverse( $copy->raw )) ) {
        die( __PACKAGE__.": error reversing vector " . gsl_strerror($status) );
    }
    return $copy;
}

=head2  set() 

Sets values of an of a Math::GSL::Vector object.

    my $vector = Math::GSL::VectorComplex->new(3);
    $vector->set([1,2], [8,23]);

This sets the second and third value to 8 and 23.

=cut

sub set {
    my ($self, $indices, $values) = @_;
    die (__PACKAGE__.'::set($indices, $values) - $indices and $values must be array references of the same length')
        unless ( ref $indices eq 'ARRAY' && ref $values eq 'ARRAY' &&  $#$indices == $#$values );
    eval {
        map {  gsl_vector_complex_set($self->{_vector}, $indices->[$_], $values->[$_] ) } (0..$#$indices);
    };
    # better error handling?
    warn $@ if $@;
    return;
}

=head2 copy()

Returns a copy of the vector, which has the same length and values but resides at a different location in memory.

    my $vector = Math::GSL::VectorComplex->new([10 .. 20]);
    my $copy   = $vector->copy;

=cut

sub copy {
    my $self = shift;
    my $copy = Math::GSL::VectorComplex->new( $self->length );
    my $status = gsl_vector_complex_memcpy($copy->raw, $self->raw);
    if ( $status != $GSL_SUCCESS ) {
        croak "Math::GSL - error copying memory, aborting. $! status=$status";
    }
    return $copy;
}

=head2 swap()

Exchanges the values in the vectors $v with $w by copying.

    my $v = Math::GSL::VectorComplex->new([1..5]);
    my $w = Math::GSL::VectorComplex->new([3..7]);
    $v->swap( $w );

=cut

sub swap() {
    my ($self,$other) = @_;
    croak "Math::GSL::VectorComplex : \$v->swap(\$w) - \$w must be a Math::GSL::VectorComplex"
        unless ref $other eq 'Math::GSL::VectorComplex';
    gsl_vector_complex_swap( $self->raw, $other->raw );
    return $self;
}

sub _multiplication {
    my ($left,$right) = @_;
    my $lcopy = $left->copy;

    if ( blessed $right && $right->isa(__PACKAGE__) ) {
        return $lcopy->dot_product($right);
    } else {
        # will be in upcoming gsl 1.12
        # gsl_vector_complex_scale($lcopy->raw, $right);
    }
    return $lcopy;
}

sub _subtract {
    my ($left, $right, $flip) = @_;

    if ($flip) {
        my $lcopy = $left->copy;
        # will be in upcoming gsl 1.12
        # gsl_vector_complex_scale($lcopy->raw, -1 );
        gsl_vector_add_constant($lcopy->raw, $right);
        return $lcopy;
    } else {
        return _addition($left, -1.0*$right);
    }
}

sub _addition {
    my ($left, $right, $flip) = @_;
    my $lcopy = $left->copy;

    if ( blessed $right && $right->isa('Math::GSL::Vector') && blessed $left && $left->isa('Math::GSL::Vector') ) {
        if ( $left->length == $right->length ) {
            gsl_vector_complex_add($lcopy->raw, $right->raw);
        } else {
            croak "Math::GSL - addition of vectors must be called with two objects vectors and must have the same length";
        }
    } else {
        gsl_vector_complex_add_constant($lcopy->raw, $right);
    }
    return $lcopy;
}

sub dot_product_pp {
    my ($left,$right) = @_;
    my $sum=0;
    if ( blessed $right && $right->isa('Math::GSL::Vector') && $left->length == $right->length ) {
         my @l = $left->as_list;
         my @r = $right->as_list;
         map { $sum += $l[$_] * $r[$_] } (0..$#l);
        return $sum;
    } else {
        croak "dot_product() must be called with two vectors";
    }
}

sub dot_product {
    my ($left,$right) = @_;

    my ($status, $product) = gsl_blas_ddot($left->raw,$right->raw);
    croak sprintf "Math::GSL::dot_product - %s", gsl_strerror($status) if ($status != $GSL_SUCCESS);
    return $product;
}

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

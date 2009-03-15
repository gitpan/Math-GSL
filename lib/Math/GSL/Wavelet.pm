# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Wavelet;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Waveletc;
bootstrap Math::GSL::Wavelet;
package Math::GSL::Wavelet;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Wavelet;

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

package Math::GSL::Wavelet;

*gsl_wavelet_alloc = *Math::GSL::Waveletc::gsl_wavelet_alloc;
*gsl_wavelet_free = *Math::GSL::Waveletc::gsl_wavelet_free;
*gsl_wavelet_name = *Math::GSL::Waveletc::gsl_wavelet_name;
*gsl_wavelet_workspace_alloc = *Math::GSL::Waveletc::gsl_wavelet_workspace_alloc;
*gsl_wavelet_workspace_free = *Math::GSL::Waveletc::gsl_wavelet_workspace_free;
*gsl_wavelet_transform = *Math::GSL::Waveletc::gsl_wavelet_transform;
*gsl_wavelet_transform_forward = *Math::GSL::Waveletc::gsl_wavelet_transform_forward;
*gsl_wavelet_transform_inverse = *Math::GSL::Waveletc::gsl_wavelet_transform_inverse;

############# Class : Math::GSL::Wavelet::gsl_wavelet_type ##############

package Math::GSL::Wavelet::gsl_wavelet_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Wavelet );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Waveletc::gsl_wavelet_type_name_get;
*swig_name_set = *Math::GSL::Waveletc::gsl_wavelet_type_name_set;
*swig_init_get = *Math::GSL::Waveletc::gsl_wavelet_type_init_get;
*swig_init_set = *Math::GSL::Waveletc::gsl_wavelet_type_init_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Waveletc::new_gsl_wavelet_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Waveletc::delete_gsl_wavelet_type($self);
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


############# Class : Math::GSL::Wavelet::gsl_wavelet ##############

package Math::GSL::Wavelet::gsl_wavelet;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Wavelet );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Waveletc::gsl_wavelet_type_get;
*swig_type_set = *Math::GSL::Waveletc::gsl_wavelet_type_set;
*swig_h1_get = *Math::GSL::Waveletc::gsl_wavelet_h1_get;
*swig_h1_set = *Math::GSL::Waveletc::gsl_wavelet_h1_set;
*swig_g1_get = *Math::GSL::Waveletc::gsl_wavelet_g1_get;
*swig_g1_set = *Math::GSL::Waveletc::gsl_wavelet_g1_set;
*swig_h2_get = *Math::GSL::Waveletc::gsl_wavelet_h2_get;
*swig_h2_set = *Math::GSL::Waveletc::gsl_wavelet_h2_set;
*swig_g2_get = *Math::GSL::Waveletc::gsl_wavelet_g2_get;
*swig_g2_set = *Math::GSL::Waveletc::gsl_wavelet_g2_set;
*swig_nc_get = *Math::GSL::Waveletc::gsl_wavelet_nc_get;
*swig_nc_set = *Math::GSL::Waveletc::gsl_wavelet_nc_set;
*swig_offset_get = *Math::GSL::Waveletc::gsl_wavelet_offset_get;
*swig_offset_set = *Math::GSL::Waveletc::gsl_wavelet_offset_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Waveletc::new_gsl_wavelet(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Waveletc::delete_gsl_wavelet($self);
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


############# Class : Math::GSL::Wavelet::gsl_wavelet_workspace ##############

package Math::GSL::Wavelet::gsl_wavelet_workspace;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Wavelet );
%OWNER = ();
%ITERATORS = ();
*swig_scratch_get = *Math::GSL::Waveletc::gsl_wavelet_workspace_scratch_get;
*swig_scratch_set = *Math::GSL::Waveletc::gsl_wavelet_workspace_scratch_set;
*swig_n_get = *Math::GSL::Waveletc::gsl_wavelet_workspace_n_get;
*swig_n_set = *Math::GSL::Waveletc::gsl_wavelet_workspace_n_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Waveletc::new_gsl_wavelet_workspace(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Waveletc::delete_gsl_wavelet_workspace($self);
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

package Math::GSL::Wavelet;

*GSL_MAJOR_VERSION = *Math::GSL::Waveletc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Waveletc::GSL_MINOR_VERSION;
*forward = *Math::GSL::Waveletc::forward;
*backward = *Math::GSL::Waveletc::backward;
*gsl_wavelet_forward = *Math::GSL::Waveletc::gsl_wavelet_forward;
*gsl_wavelet_backward = *Math::GSL::Waveletc::gsl_wavelet_backward;

my %__gsl_wavelet_daubechies_hash;
tie %__gsl_wavelet_daubechies_hash,"Math::GSL::Wavelet::gsl_wavelet_type", $Math::GSL::Waveletc::gsl_wavelet_daubechies;
$gsl_wavelet_daubechies= \%__gsl_wavelet_daubechies_hash;
bless $gsl_wavelet_daubechies, Math::GSL::Wavelet::gsl_wavelet_type;

my %__gsl_wavelet_daubechies_centered_hash;
tie %__gsl_wavelet_daubechies_centered_hash,"Math::GSL::Wavelet::gsl_wavelet_type", $Math::GSL::Waveletc::gsl_wavelet_daubechies_centered;
$gsl_wavelet_daubechies_centered= \%__gsl_wavelet_daubechies_centered_hash;
bless $gsl_wavelet_daubechies_centered, Math::GSL::Wavelet::gsl_wavelet_type;

my %__gsl_wavelet_haar_hash;
tie %__gsl_wavelet_haar_hash,"Math::GSL::Wavelet::gsl_wavelet_type", $Math::GSL::Waveletc::gsl_wavelet_haar;
$gsl_wavelet_haar= \%__gsl_wavelet_haar_hash;
bless $gsl_wavelet_haar, Math::GSL::Wavelet::gsl_wavelet_type;

my %__gsl_wavelet_haar_centered_hash;
tie %__gsl_wavelet_haar_centered_hash,"Math::GSL::Wavelet::gsl_wavelet_type", $Math::GSL::Waveletc::gsl_wavelet_haar_centered;
$gsl_wavelet_haar_centered= \%__gsl_wavelet_haar_centered_hash;
bless $gsl_wavelet_haar_centered, Math::GSL::Wavelet::gsl_wavelet_type;

my %__gsl_wavelet_bspline_hash;
tie %__gsl_wavelet_bspline_hash,"Math::GSL::Wavelet::gsl_wavelet_type", $Math::GSL::Waveletc::gsl_wavelet_bspline;
$gsl_wavelet_bspline= \%__gsl_wavelet_bspline_hash;
bless $gsl_wavelet_bspline, Math::GSL::Wavelet::gsl_wavelet_type;

my %__gsl_wavelet_bspline_centered_hash;
tie %__gsl_wavelet_bspline_centered_hash,"Math::GSL::Wavelet::gsl_wavelet_type", $Math::GSL::Waveletc::gsl_wavelet_bspline_centered;
$gsl_wavelet_bspline_centered= \%__gsl_wavelet_bspline_centered_hash;
bless $gsl_wavelet_bspline_centered, Math::GSL::Wavelet::gsl_wavelet_type;



@EXPORT_OK = qw/
    gsl_wavelet_alloc 
    gsl_wavelet_free 
    gsl_wavelet_name 
    gsl_wavelet_workspace_alloc 
    gsl_wavelet_workspace_free 
    gsl_wavelet_transform 
    gsl_wavelet_transform_forward 
    gsl_wavelet_transform_inverse 
    $gsl_wavelet_daubechies
    $gsl_wavelet_daubechies_centered
    $gsl_wavelet_haar
    $gsl_wavelet_haar_centered
    $gsl_wavelet_bspline
    $gsl_wavelet_bspline_centered
/;


%EXPORT_TAGS = ( all => [ @EXPORT_OK ], );

=head1 NAME

Math::GSL::Wavelet - Wavelets (for 1-D real data)

=head1 SYNOPSIS

    use Math::GSL::Wavelet qw/:all/;

=cut

=head1 DESCRIPTION

Here is a list of all the functions included in this module :

=over 1

=item C<gsl_wavelet_alloc($T, $k)> - This function allocates and initializes a wavelet object of type $T, where $T must be one of the constants below. The parameter $k selects the specific member of the wavelet family.

=item C<gsl_wavelet_free($w)> - This function frees the wavelet object $w. 

=item C<gsl_wavelet_name>

=item C<gsl_wavelet_workspace_alloc($n)> - This function allocates a workspace for the discrete wavelet transform. To perform a one-dimensional transform on $n elements, a workspace of size $n must be provided. For two-dimensional transforms of $n-by-$n matrices it is sufficient to allocate a workspace of size $n, since the transform operates on individual rows and columns.

=item C<gsl_wavelet_workspace_free($work)> - This function frees the allocated workspace work. 

=item C<gsl_wavelet_transform>

=item C<gsl_wavelet_transform_forward($w, $data, $stride, $n, $work)> - This functions compute in-place forward discrete wavelet transforms of length $n with stride $stride on the array $data. The length of the transform $n is restricted to powers of two. For the forward transform, the elements of the original array are replaced by the discrete wavelet transform f_i -> w_{j,k} in a packed triangular storage layout, where j is the index of the level j = 0 ... J-1 and k is the index of the coefficient within each level, k = 0 ... (2^j)-1. The total number of levels is J = \log_2(n). The output data has the following form,

=over

=item (s_{-1,0}, d_{0,0}, d_{1,0}, d_{1,1}, d_{2,0}, ...,

=item d_{j,k}, ..., d_{J-1,2^{J-1}-1})

=back
     
where the first element is the smoothing coefficient s_{-1,0}, followed by the detail coefficients d_{j,k} for each level j. The backward transform inverts these coefficients to obtain the original data. These functions return a status of $GSL_SUCCESS upon successful completion. $GSL_EINVAL is returned if $n is not an integer power of 2 or if insufficient workspace is provided.

=item C<gsl_wavelet_transform_inverse>

=back

This module also contains the following constants with their valid k value for the gsl_wavelet_alloc function :

=over 1

=item $gsl_wavelet_daubechies

=item $gsl_wavelet_daubechies_centered

=back

This is the Daubechies wavelet family of maximum phase with k/2 vanishing moments. The implemented wavelets are k=4, 6, ..., 20, with k even.

=over 1

=item $gsl_wavelet_haar

=item $gsl_wavelet_haar_centered

=back

This is the Haar wavelet. The only valid choice of k for the Haar wavelet is k=2. 

=over 1

=item $gsl_wavelet_bspline

=item $gsl_wavelet_bspline_centered

=back

This is the biorthogonal B-spline wavelet family of order (i,j). The implemented values of k = 100*i + j are 103, 105, 202, 204, 206, 208, 301, 303, 305 307, 309.

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
1;

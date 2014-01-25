# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::FFT;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::FFTc;
bootstrap Math::GSL::FFT;
package Math::GSL::FFT;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::FFT;

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

package Math::GSL::FFT;

*gsl_log1p = *Math::GSL::FFTc::gsl_log1p;
*gsl_expm1 = *Math::GSL::FFTc::gsl_expm1;
*gsl_hypot = *Math::GSL::FFTc::gsl_hypot;
*gsl_hypot3 = *Math::GSL::FFTc::gsl_hypot3;
*gsl_acosh = *Math::GSL::FFTc::gsl_acosh;
*gsl_asinh = *Math::GSL::FFTc::gsl_asinh;
*gsl_atanh = *Math::GSL::FFTc::gsl_atanh;
*gsl_isnan = *Math::GSL::FFTc::gsl_isnan;
*gsl_isinf = *Math::GSL::FFTc::gsl_isinf;
*gsl_finite = *Math::GSL::FFTc::gsl_finite;
*gsl_nan = *Math::GSL::FFTc::gsl_nan;
*gsl_posinf = *Math::GSL::FFTc::gsl_posinf;
*gsl_neginf = *Math::GSL::FFTc::gsl_neginf;
*gsl_fdiv = *Math::GSL::FFTc::gsl_fdiv;
*gsl_coerce_double = *Math::GSL::FFTc::gsl_coerce_double;
*gsl_coerce_float = *Math::GSL::FFTc::gsl_coerce_float;
*gsl_coerce_long_double = *Math::GSL::FFTc::gsl_coerce_long_double;
*gsl_ldexp = *Math::GSL::FFTc::gsl_ldexp;
*gsl_frexp = *Math::GSL::FFTc::gsl_frexp;
*gsl_fcmp = *Math::GSL::FFTc::gsl_fcmp;
*gsl_pow_2 = *Math::GSL::FFTc::gsl_pow_2;
*gsl_pow_3 = *Math::GSL::FFTc::gsl_pow_3;
*gsl_pow_4 = *Math::GSL::FFTc::gsl_pow_4;
*gsl_pow_5 = *Math::GSL::FFTc::gsl_pow_5;
*gsl_pow_6 = *Math::GSL::FFTc::gsl_pow_6;
*gsl_pow_7 = *Math::GSL::FFTc::gsl_pow_7;
*gsl_pow_8 = *Math::GSL::FFTc::gsl_pow_8;
*gsl_pow_9 = *Math::GSL::FFTc::gsl_pow_9;
*gsl_pow_int = *Math::GSL::FFTc::gsl_pow_int;
*gsl_pow_uint = *Math::GSL::FFTc::gsl_pow_uint;
*gsl_fft_complex_radix2_forward = *Math::GSL::FFTc::gsl_fft_complex_radix2_forward;
*gsl_fft_complex_radix2_backward = *Math::GSL::FFTc::gsl_fft_complex_radix2_backward;
*gsl_fft_complex_radix2_inverse = *Math::GSL::FFTc::gsl_fft_complex_radix2_inverse;
*gsl_fft_complex_radix2_transform = *Math::GSL::FFTc::gsl_fft_complex_radix2_transform;
*gsl_fft_complex_radix2_dif_forward = *Math::GSL::FFTc::gsl_fft_complex_radix2_dif_forward;
*gsl_fft_complex_radix2_dif_backward = *Math::GSL::FFTc::gsl_fft_complex_radix2_dif_backward;
*gsl_fft_complex_radix2_dif_inverse = *Math::GSL::FFTc::gsl_fft_complex_radix2_dif_inverse;
*gsl_fft_complex_radix2_dif_transform = *Math::GSL::FFTc::gsl_fft_complex_radix2_dif_transform;
*gsl_fft_complex_wavetable_alloc = *Math::GSL::FFTc::gsl_fft_complex_wavetable_alloc;
*gsl_fft_complex_wavetable_free = *Math::GSL::FFTc::gsl_fft_complex_wavetable_free;
*gsl_fft_complex_workspace_alloc = *Math::GSL::FFTc::gsl_fft_complex_workspace_alloc;
*gsl_fft_complex_workspace_free = *Math::GSL::FFTc::gsl_fft_complex_workspace_free;
*gsl_fft_complex_memcpy = *Math::GSL::FFTc::gsl_fft_complex_memcpy;
*gsl_fft_complex_forward = *Math::GSL::FFTc::gsl_fft_complex_forward;
*gsl_fft_complex_backward = *Math::GSL::FFTc::gsl_fft_complex_backward;
*gsl_fft_complex_inverse = *Math::GSL::FFTc::gsl_fft_complex_inverse;
*gsl_fft_complex_transform = *Math::GSL::FFTc::gsl_fft_complex_transform;
*gsl_fft_halfcomplex_radix2_backward = *Math::GSL::FFTc::gsl_fft_halfcomplex_radix2_backward;
*gsl_fft_halfcomplex_radix2_inverse = *Math::GSL::FFTc::gsl_fft_halfcomplex_radix2_inverse;
*gsl_fft_halfcomplex_radix2_transform = *Math::GSL::FFTc::gsl_fft_halfcomplex_radix2_transform;
*gsl_fft_halfcomplex_wavetable_alloc = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_alloc;
*gsl_fft_halfcomplex_wavetable_free = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_free;
*gsl_fft_halfcomplex_backward = *Math::GSL::FFTc::gsl_fft_halfcomplex_backward;
*gsl_fft_halfcomplex_inverse = *Math::GSL::FFTc::gsl_fft_halfcomplex_inverse;
*gsl_fft_halfcomplex_transform = *Math::GSL::FFTc::gsl_fft_halfcomplex_transform;
*gsl_fft_halfcomplex_unpack = *Math::GSL::FFTc::gsl_fft_halfcomplex_unpack;
*gsl_fft_halfcomplex_radix2_unpack = *Math::GSL::FFTc::gsl_fft_halfcomplex_radix2_unpack;
*gsl_fft_real_radix2_transform = *Math::GSL::FFTc::gsl_fft_real_radix2_transform;
*gsl_fft_real_wavetable_alloc = *Math::GSL::FFTc::gsl_fft_real_wavetable_alloc;
*gsl_fft_real_wavetable_free = *Math::GSL::FFTc::gsl_fft_real_wavetable_free;
*gsl_fft_real_workspace_alloc = *Math::GSL::FFTc::gsl_fft_real_workspace_alloc;
*gsl_fft_real_workspace_free = *Math::GSL::FFTc::gsl_fft_real_workspace_free;
*gsl_fft_real_transform = *Math::GSL::FFTc::gsl_fft_real_transform;
*gsl_fft_real_unpack = *Math::GSL::FFTc::gsl_fft_real_unpack;

############# Class : Math::GSL::FFT::gsl_function_struct ##############

package Math::GSL::FFT::gsl_function_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::FFTc::gsl_function_struct_function_get;
*swig_function_set = *Math::GSL::FFTc::gsl_function_struct_function_set;
*swig_params_get = *Math::GSL::FFTc::gsl_function_struct_params_get;
*swig_params_set = *Math::GSL::FFTc::gsl_function_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_function_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_function_struct($self);
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


############# Class : Math::GSL::FFT::gsl_function_fdf_struct ##############

package Math::GSL::FFT::gsl_function_fdf_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::FFTc::gsl_function_fdf_struct_f_get;
*swig_f_set = *Math::GSL::FFTc::gsl_function_fdf_struct_f_set;
*swig_df_get = *Math::GSL::FFTc::gsl_function_fdf_struct_df_get;
*swig_df_set = *Math::GSL::FFTc::gsl_function_fdf_struct_df_set;
*swig_fdf_get = *Math::GSL::FFTc::gsl_function_fdf_struct_fdf_get;
*swig_fdf_set = *Math::GSL::FFTc::gsl_function_fdf_struct_fdf_set;
*swig_params_get = *Math::GSL::FFTc::gsl_function_fdf_struct_params_get;
*swig_params_set = *Math::GSL::FFTc::gsl_function_fdf_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_function_fdf_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_function_fdf_struct($self);
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


############# Class : Math::GSL::FFT::gsl_function_vec_struct ##############

package Math::GSL::FFT::gsl_function_vec_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::FFTc::gsl_function_vec_struct_function_get;
*swig_function_set = *Math::GSL::FFTc::gsl_function_vec_struct_function_set;
*swig_params_get = *Math::GSL::FFTc::gsl_function_vec_struct_params_get;
*swig_params_set = *Math::GSL::FFTc::gsl_function_vec_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_function_vec_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_function_vec_struct($self);
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


############# Class : Math::GSL::FFT::gsl_complex_long_double ##############

package Math::GSL::FFT::gsl_complex_long_double;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::FFTc::gsl_complex_long_double_dat_get;
*swig_dat_set = *Math::GSL::FFTc::gsl_complex_long_double_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_complex_long_double(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_complex_long_double($self);
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


############# Class : Math::GSL::FFT::gsl_complex ##############

package Math::GSL::FFT::gsl_complex;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::FFTc::gsl_complex_dat_get;
*swig_dat_set = *Math::GSL::FFTc::gsl_complex_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_complex(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_complex($self);
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


############# Class : Math::GSL::FFT::gsl_complex_float ##############

package Math::GSL::FFT::gsl_complex_float;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_dat_get = *Math::GSL::FFTc::gsl_complex_float_dat_get;
*swig_dat_set = *Math::GSL::FFTc::gsl_complex_float_dat_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_complex_float(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_complex_float($self);
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


############# Class : Math::GSL::FFT::gsl_fft_complex_wavetable ##############

package Math::GSL::FFT::gsl_fft_complex_wavetable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::FFTc::gsl_fft_complex_wavetable_n_get;
*swig_n_set = *Math::GSL::FFTc::gsl_fft_complex_wavetable_n_set;
*swig_nf_get = *Math::GSL::FFTc::gsl_fft_complex_wavetable_nf_get;
*swig_nf_set = *Math::GSL::FFTc::gsl_fft_complex_wavetable_nf_set;
*swig_factor_get = *Math::GSL::FFTc::gsl_fft_complex_wavetable_factor_get;
*swig_factor_set = *Math::GSL::FFTc::gsl_fft_complex_wavetable_factor_set;
*swig_twiddle_get = *Math::GSL::FFTc::gsl_fft_complex_wavetable_twiddle_get;
*swig_twiddle_set = *Math::GSL::FFTc::gsl_fft_complex_wavetable_twiddle_set;
*swig_trig_get = *Math::GSL::FFTc::gsl_fft_complex_wavetable_trig_get;
*swig_trig_set = *Math::GSL::FFTc::gsl_fft_complex_wavetable_trig_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_fft_complex_wavetable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_fft_complex_wavetable($self);
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


############# Class : Math::GSL::FFT::gsl_fft_complex_workspace ##############

package Math::GSL::FFT::gsl_fft_complex_workspace;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::FFTc::gsl_fft_complex_workspace_n_get;
*swig_n_set = *Math::GSL::FFTc::gsl_fft_complex_workspace_n_set;
*swig_scratch_get = *Math::GSL::FFTc::gsl_fft_complex_workspace_scratch_get;
*swig_scratch_set = *Math::GSL::FFTc::gsl_fft_complex_workspace_scratch_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_fft_complex_workspace(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_fft_complex_workspace($self);
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


############# Class : Math::GSL::FFT::gsl_fft_halfcomplex_wavetable ##############

package Math::GSL::FFT::gsl_fft_halfcomplex_wavetable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_n_get;
*swig_n_set = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_n_set;
*swig_nf_get = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_nf_get;
*swig_nf_set = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_nf_set;
*swig_factor_get = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_factor_get;
*swig_factor_set = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_factor_set;
*swig_twiddle_get = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_twiddle_get;
*swig_twiddle_set = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_twiddle_set;
*swig_trig_get = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_trig_get;
*swig_trig_set = *Math::GSL::FFTc::gsl_fft_halfcomplex_wavetable_trig_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_fft_halfcomplex_wavetable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_fft_halfcomplex_wavetable($self);
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


############# Class : Math::GSL::FFT::gsl_fft_real_wavetable ##############

package Math::GSL::FFT::gsl_fft_real_wavetable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::FFTc::gsl_fft_real_wavetable_n_get;
*swig_n_set = *Math::GSL::FFTc::gsl_fft_real_wavetable_n_set;
*swig_nf_get = *Math::GSL::FFTc::gsl_fft_real_wavetable_nf_get;
*swig_nf_set = *Math::GSL::FFTc::gsl_fft_real_wavetable_nf_set;
*swig_factor_get = *Math::GSL::FFTc::gsl_fft_real_wavetable_factor_get;
*swig_factor_set = *Math::GSL::FFTc::gsl_fft_real_wavetable_factor_set;
*swig_twiddle_get = *Math::GSL::FFTc::gsl_fft_real_wavetable_twiddle_get;
*swig_twiddle_set = *Math::GSL::FFTc::gsl_fft_real_wavetable_twiddle_set;
*swig_trig_get = *Math::GSL::FFTc::gsl_fft_real_wavetable_trig_get;
*swig_trig_set = *Math::GSL::FFTc::gsl_fft_real_wavetable_trig_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_fft_real_wavetable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_fft_real_wavetable($self);
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


############# Class : Math::GSL::FFT::gsl_fft_real_workspace ##############

package Math::GSL::FFT::gsl_fft_real_workspace;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::FFT );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::FFTc::gsl_fft_real_workspace_n_get;
*swig_n_set = *Math::GSL::FFTc::gsl_fft_real_workspace_n_set;
*swig_scratch_get = *Math::GSL::FFTc::gsl_fft_real_workspace_scratch_get;
*swig_scratch_set = *Math::GSL::FFTc::gsl_fft_real_workspace_scratch_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::FFTc::new_gsl_fft_real_workspace(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::FFTc::delete_gsl_fft_real_workspace($self);
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

package Math::GSL::FFT;

*GSL_MAJOR_VERSION = *Math::GSL::FFTc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::FFTc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::FFTc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::FFTc::GSL_NEGZERO;
*M_E = *Math::GSL::FFTc::M_E;
*M_LOG2E = *Math::GSL::FFTc::M_LOG2E;
*M_LOG10E = *Math::GSL::FFTc::M_LOG10E;
*M_SQRT2 = *Math::GSL::FFTc::M_SQRT2;
*M_SQRT1_2 = *Math::GSL::FFTc::M_SQRT1_2;
*M_SQRT3 = *Math::GSL::FFTc::M_SQRT3;
*M_PI = *Math::GSL::FFTc::M_PI;
*M_PI_2 = *Math::GSL::FFTc::M_PI_2;
*M_PI_4 = *Math::GSL::FFTc::M_PI_4;
*M_SQRTPI = *Math::GSL::FFTc::M_SQRTPI;
*M_2_SQRTPI = *Math::GSL::FFTc::M_2_SQRTPI;
*M_1_PI = *Math::GSL::FFTc::M_1_PI;
*M_2_PI = *Math::GSL::FFTc::M_2_PI;
*M_LN10 = *Math::GSL::FFTc::M_LN10;
*M_LN2 = *Math::GSL::FFTc::M_LN2;
*M_LNPI = *Math::GSL::FFTc::M_LNPI;
*M_EULER = *Math::GSL::FFTc::M_EULER;
*GSL_DBL_EPSILON = *Math::GSL::FFTc::GSL_DBL_EPSILON;
*GSL_SQRT_DBL_EPSILON = *Math::GSL::FFTc::GSL_SQRT_DBL_EPSILON;
*GSL_ROOT3_DBL_EPSILON = *Math::GSL::FFTc::GSL_ROOT3_DBL_EPSILON;
*GSL_ROOT4_DBL_EPSILON = *Math::GSL::FFTc::GSL_ROOT4_DBL_EPSILON;
*GSL_ROOT5_DBL_EPSILON = *Math::GSL::FFTc::GSL_ROOT5_DBL_EPSILON;
*GSL_ROOT6_DBL_EPSILON = *Math::GSL::FFTc::GSL_ROOT6_DBL_EPSILON;
*GSL_LOG_DBL_EPSILON = *Math::GSL::FFTc::GSL_LOG_DBL_EPSILON;
*GSL_DBL_MIN = *Math::GSL::FFTc::GSL_DBL_MIN;
*GSL_SQRT_DBL_MIN = *Math::GSL::FFTc::GSL_SQRT_DBL_MIN;
*GSL_ROOT3_DBL_MIN = *Math::GSL::FFTc::GSL_ROOT3_DBL_MIN;
*GSL_ROOT4_DBL_MIN = *Math::GSL::FFTc::GSL_ROOT4_DBL_MIN;
*GSL_ROOT5_DBL_MIN = *Math::GSL::FFTc::GSL_ROOT5_DBL_MIN;
*GSL_ROOT6_DBL_MIN = *Math::GSL::FFTc::GSL_ROOT6_DBL_MIN;
*GSL_LOG_DBL_MIN = *Math::GSL::FFTc::GSL_LOG_DBL_MIN;
*GSL_DBL_MAX = *Math::GSL::FFTc::GSL_DBL_MAX;
*GSL_SQRT_DBL_MAX = *Math::GSL::FFTc::GSL_SQRT_DBL_MAX;
*GSL_ROOT3_DBL_MAX = *Math::GSL::FFTc::GSL_ROOT3_DBL_MAX;
*GSL_ROOT4_DBL_MAX = *Math::GSL::FFTc::GSL_ROOT4_DBL_MAX;
*GSL_ROOT5_DBL_MAX = *Math::GSL::FFTc::GSL_ROOT5_DBL_MAX;
*GSL_ROOT6_DBL_MAX = *Math::GSL::FFTc::GSL_ROOT6_DBL_MAX;
*GSL_LOG_DBL_MAX = *Math::GSL::FFTc::GSL_LOG_DBL_MAX;
*GSL_FLT_EPSILON = *Math::GSL::FFTc::GSL_FLT_EPSILON;
*GSL_SQRT_FLT_EPSILON = *Math::GSL::FFTc::GSL_SQRT_FLT_EPSILON;
*GSL_ROOT3_FLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT3_FLT_EPSILON;
*GSL_ROOT4_FLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT4_FLT_EPSILON;
*GSL_ROOT5_FLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT5_FLT_EPSILON;
*GSL_ROOT6_FLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT6_FLT_EPSILON;
*GSL_LOG_FLT_EPSILON = *Math::GSL::FFTc::GSL_LOG_FLT_EPSILON;
*GSL_FLT_MIN = *Math::GSL::FFTc::GSL_FLT_MIN;
*GSL_SQRT_FLT_MIN = *Math::GSL::FFTc::GSL_SQRT_FLT_MIN;
*GSL_ROOT3_FLT_MIN = *Math::GSL::FFTc::GSL_ROOT3_FLT_MIN;
*GSL_ROOT4_FLT_MIN = *Math::GSL::FFTc::GSL_ROOT4_FLT_MIN;
*GSL_ROOT5_FLT_MIN = *Math::GSL::FFTc::GSL_ROOT5_FLT_MIN;
*GSL_ROOT6_FLT_MIN = *Math::GSL::FFTc::GSL_ROOT6_FLT_MIN;
*GSL_LOG_FLT_MIN = *Math::GSL::FFTc::GSL_LOG_FLT_MIN;
*GSL_FLT_MAX = *Math::GSL::FFTc::GSL_FLT_MAX;
*GSL_SQRT_FLT_MAX = *Math::GSL::FFTc::GSL_SQRT_FLT_MAX;
*GSL_ROOT3_FLT_MAX = *Math::GSL::FFTc::GSL_ROOT3_FLT_MAX;
*GSL_ROOT4_FLT_MAX = *Math::GSL::FFTc::GSL_ROOT4_FLT_MAX;
*GSL_ROOT5_FLT_MAX = *Math::GSL::FFTc::GSL_ROOT5_FLT_MAX;
*GSL_ROOT6_FLT_MAX = *Math::GSL::FFTc::GSL_ROOT6_FLT_MAX;
*GSL_LOG_FLT_MAX = *Math::GSL::FFTc::GSL_LOG_FLT_MAX;
*GSL_SFLT_EPSILON = *Math::GSL::FFTc::GSL_SFLT_EPSILON;
*GSL_SQRT_SFLT_EPSILON = *Math::GSL::FFTc::GSL_SQRT_SFLT_EPSILON;
*GSL_ROOT3_SFLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT3_SFLT_EPSILON;
*GSL_ROOT4_SFLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT4_SFLT_EPSILON;
*GSL_ROOT5_SFLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT5_SFLT_EPSILON;
*GSL_ROOT6_SFLT_EPSILON = *Math::GSL::FFTc::GSL_ROOT6_SFLT_EPSILON;
*GSL_LOG_SFLT_EPSILON = *Math::GSL::FFTc::GSL_LOG_SFLT_EPSILON;
*GSL_MACH_EPS = *Math::GSL::FFTc::GSL_MACH_EPS;
*GSL_SQRT_MACH_EPS = *Math::GSL::FFTc::GSL_SQRT_MACH_EPS;
*GSL_ROOT3_MACH_EPS = *Math::GSL::FFTc::GSL_ROOT3_MACH_EPS;
*GSL_ROOT4_MACH_EPS = *Math::GSL::FFTc::GSL_ROOT4_MACH_EPS;
*GSL_ROOT5_MACH_EPS = *Math::GSL::FFTc::GSL_ROOT5_MACH_EPS;
*GSL_ROOT6_MACH_EPS = *Math::GSL::FFTc::GSL_ROOT6_MACH_EPS;
*GSL_LOG_MACH_EPS = *Math::GSL::FFTc::GSL_LOG_MACH_EPS;
*gsl_fft_forward = *Math::GSL::FFTc::gsl_fft_forward;
*gsl_fft_backward = *Math::GSL::FFTc::gsl_fft_backward;

@EXPORT_complex = qw/
               gsl_fft_complex_radix2_forward 
               gsl_fft_complex_radix2_backward 
               gsl_fft_complex_radix2_inverse 
               gsl_fft_complex_radix2_transform 
               gsl_fft_complex_radix2_dif_forward 
               gsl_fft_complex_radix2_dif_backward 
               gsl_fft_complex_radix2_dif_inverse 
               gsl_fft_complex_radix2_dif_transform 
               gsl_fft_complex_wavetable_alloc 
               gsl_fft_complex_wavetable_free 
               gsl_fft_complex_workspace_alloc 
               gsl_fft_complex_workspace_free 
               gsl_fft_complex_memcpy 
               gsl_fft_complex_forward 
               gsl_fft_complex_backward 
               gsl_fft_complex_inverse 
               gsl_fft_complex_transform 
               /;
@EXPORT_halfcomplex = qw/
               gsl_fft_halfcomplex_radix2_backward 
               gsl_fft_halfcomplex_radix2_inverse 
               gsl_fft_halfcomplex_radix2_transform 
               gsl_fft_halfcomplex_wavetable_alloc 
               gsl_fft_halfcomplex_wavetable_free 
               gsl_fft_halfcomplex_backward 
               gsl_fft_halfcomplex_inverse 
               gsl_fft_halfcomplex_transform 
               gsl_fft_halfcomplex_unpack 
               gsl_fft_halfcomplex_radix2_unpack 
               /;
@EXPORT_real = qw/ 
               gsl_fft_real_radix2_transform 
               gsl_fft_real_wavetable_alloc 
               gsl_fft_real_wavetable_free 
               gsl_fft_real_workspace_alloc 
               gsl_fft_real_workspace_free 
               gsl_fft_real_transform 
               gsl_fft_real_unpack 
             /;
@EXPORT_vars = qw/
                $gsl_fft_forward
                $gsl_fft_backward
                /;
@EXPORT_OK =   (
                @EXPORT_real, 
                @EXPORT_complex,
                @EXPORT_halfcomplex, 
                @EXPORT_vars,
                );
%EXPORT_TAGS = ( 
                all         => \@EXPORT_OK, 
                real        => \@EXPORT_real,
                complex     => \@EXPORT_complex,
                halfcomplex => \@EXPORT_halfcomplex,
                vars        => \@EXPORT_vars,
               );
__END__

=head1 NAME

Math::GSL::FFT - Fast Fourier Transforms (FFT)

=head1 SYNOPSIS

    use Math::GSL::FFT qw /:all/;
    # alternating elements are real/imaginary part, hence 256 element array
    my $data = [ (1) x 10, (0) x 236, (1) x 10 ]; 

    # use every element of the array
    my $stride = 1;  

    # But it contains 128 complex numbers
    my ($status, $fft) = gsl_fft_complex_radix2_forward ($data, $stride, 128); 

=head1 DESCRIPTION

This module and this documentation is still in a very early state. Danger Will Robinson!
An OO interface will evolve soon.

=over

=item * C<gsl_fft_complex_radix2_forward($data, $stride, $n) > 

This function computes the forward FFTs of length $n with stride $stride, on
the array reference $data using an in-place radix-2 decimation-in-time
algorithm. The length of the transform $n is restricted to powers of two. For
the transform version of the function the sign argument can be either forward
(-1) or backward (+1). The functions return a value of $GSL_SUCCESS if no
errors were detected, or $GSL_EDOM if the length of the data $n is not a power
of two. The complex functions of the FFT module are not yet fully implemented. 

=item * C<gsl_fft_complex_radix2_backward >

=item * C<gsl_fft_complex_radix2_inverse >

=item * C<gsl_fft_complex_radix2_transform >

=item * C<gsl_fft_complex_radix2_dif_forward >

=item * C<gsl_fft_complex_radix2_dif_backward >

=item * C<gsl_fft_complex_radix2_dif_inverse >

=item * C<gsl_fft_complex_radix2_dif_transform >

=item * C<gsl_fft_complex_wavetable_alloc($n)> 

This function prepares a trigonometric lookup table for a complex FFT of length
$n. The function returns a pointer to the newly allocated
gsl_fft_complex_wavetable if no errors were detected, and a null pointer in the
case of error. The length $n is factorized into a product of subtransforms, and
the factors and their trigonometric coefficients are stored in the wavetable.
The trigonometric coefficients are computed using direct calls to sin and cos,
for accuracy. Recursion relations could be used to compute the lookup table
faster, but if an application performs many FFTs of the same length then this
computation is a one-off overhead which does not affect the final throughput.
The wavetable structure can be used repeatedly for any transform of the same
length. The table is not modified by calls to any of the other FFT functions.
The same wavetable can be used for both forward and backward (or inverse)
transforms of a given length. 

=item * C<gsl_fft_complex_wavetable_free($wavetable)> 

This function frees the memory associated with the wavetable $wavetable. The
wavetable can be freed if no further FFTs of the same length will be needed.

=item * C<gsl_fft_complex_workspace_alloc($n)> 

This function allocates a workspace for a complex transform of length $n. 

=item * C<gsl_fft_complex_workspace_free($workspace) > 

This function frees the memory associated with the workspace $workspace. The
workspace can be freed if no further FFTs of the same length will be needed.

=item * C<gsl_fft_complex_memcpy >

=item * C<gsl_fft_complex_forward >

=item * C<gsl_fft_complex_backward >

=item * C<gsl_fft_complex_inverse >

=item * C<gsl_fft_complex_transform >

=item * C<gsl_fft_halfcomplex_radix2_backward($data, $stride, $n)> 

This function computes the backwards in-place radix-2 FFT of length $n and
stride $stride on the half-complex sequence data stored according the output
scheme used by gsl_fft_real_radix2. The result is a real array stored in
natural order.

=item * C<gsl_fft_halfcomplex_radix2_inverse($data, $stride, $n)> 

This function computes the inverse in-place radix-2 FFT of length $n and stride
$stride on the half-complex sequence data stored according the output scheme
used by gsl_fft_real_radix2. The result is a real array stored in natural
order.

=item * C<gsl_fft_halfcomplex_radix2_transform>

=item * C<gsl_fft_halfcomplex_wavetable_alloc($n)> 

This function prepares trigonometric lookup tables for an FFT of size $n real
elements. The functions return a pointer to the newly allocated struct if no
errors were detected, and a null pointer in the case of error. The length $n is
factorized into a product of subtransforms, and the factors and their
trigonometric coefficients are stored in the wavetable. The trigonometric
coefficients are computed using direct calls to sin and cos, for accuracy.
Recursion relations could be used to compute the lookup table faster, but if an
application performs many FFTs of the same length then computing the wavetable
is a one-off overhead which does not affect the final throughput.  The
wavetable structure can be used repeatedly for any transform of the same
length. The table is not modified by calls to any of the other FFT functions.
The appropriate type of wavetable must be used for forward real or inverse
half-complex transforms. 

=item * C<gsl_fft_halfcomplex_wavetable_free($wavetable)> 

This function frees the memory associated with the wavetable $wavetable. The
wavetable can be freed if no further FFTs of the same length will be needed.

=item * C<gsl_fft_halfcomplex_backward >

=item * C<gsl_fft_halfcomplex_inverse >

=item * C<gsl_fft_halfcomplex_transform >

=item * C<gsl_fft_halfcomplex_unpack >

=item * C<gsl_fft_halfcomplex_radix2_unpack >

=item * C<gsl_fft_real_radix2_transform($data, $stride, $n) > 

This function computes an in-place radix-2 FFT of length $n and stride $stride
on the real array reference $data. The output is a half-complex sequence, which
is stored in-place. The arrangement of the half-complex terms uses the
following scheme: for k < N/2 the real part of the k-th term is stored in
location k, and the corresponding imaginary part is stored in location N-k.
Terms with k > N/2 can be reconstructed using the symmetry z_k = z^*_{N-k}. The
terms for k=0 and k=N/2 are both purely real, and count as a special case.
Their real parts are stored in locations 0 and N/2 respectively, while their
imaginary parts which are zero are not stored. The following table shows the
correspondence between the output data and the equivalent results obtained by
considering the input data as a complex sequence with zero imaginary part,

          complex[0].real    =    data[0]
          complex[0].imag    =    0
          complex[1].real    =    data[1]
          complex[1].imag    =    data[N-1]
          ...............         ................
          complex[k].real    =    data[k]
          complex[k].imag    =    data[N-k]
          ...............         ................
          complex[N/2].real  =    data[N/2]
          complex[N/2].imag  =    0
          ...............         ................
          complex[k'].real   =    data[k]        k' = N - k
          complex[k'].imag   =   -data[N-k]
          ...............         ................
          complex[N-1].real  =    data[1]
          complex[N-1].imag  =   -data[N-1]

=for notyou #' 

Note that the output data can be converted into the full complex sequence using
the function gsl_fft_halfcomplex_unpack.

=item * C<gsl_fft_real_wavetable_alloc($n)> 

This function prepares trigonometric lookup tables for an FFT of size $n real
elements. The functions return a pointer to the newly allocated struct if no
errors were detected, and a null pointer in the case of error. The length $n is
factorized into a product of subtransforms, and the factors and their
trigonometric coefficients are stored in the wavetable. The trigonometric
coefficients are computed using direct calls to sin and cos, for accuracy.
Recursion relations could be used to compute the lookup table faster, but if an
application performs many FFTs of the same length then computing the wavetable
is a one-off overhead which does not affect the final throughput.  The
wavetable structure can be used repeatedly for any transform of the same
length. The table is not modified by calls to any of the other FFT functions.
The appropriate type of wavetable must be used for forward real or inverse
half-complex transforms. 

=item * C<gsl_fft_real_wavetable_free($wavetable)> 

This function frees the memory associated with the wavetable $wavetable. The
wavetable can be freed if no further FFTs of the same length will be needed. 

=item * C<gsl_fft_real_workspace_alloc($n)> 

This function allocates a workspace for a real transform of length $n. The same
workspace can be used for both forward real and inverse halfcomplex transforms.

=item * C<gsl_fft_real_workspace_free($workspace)> 

This function frees the memory associated with the workspace $workspace. The
workspace can be freed if no further FFTs of the same length will be needed.

=item * C<gsl_fft_real_transform >

=item * C<gsl_fft_real_unpack >

=back

This module also includes the following constants :

=over

=item * C<$gsl_fft_forward>

=item * C<$gsl_fft_backward>

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

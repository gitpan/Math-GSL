# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.40
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Monte;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Montec;
bootstrap Math::GSL::Monte;
package Math::GSL::Monte;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Monte;

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

package Math::GSL::Monte;

*gsl_monte_miser_integrate = *Math::GSL::Montec::gsl_monte_miser_integrate;
*gsl_monte_miser_alloc = *Math::GSL::Montec::gsl_monte_miser_alloc;
*gsl_monte_miser_init = *Math::GSL::Montec::gsl_monte_miser_init;
*gsl_monte_miser_free = *Math::GSL::Montec::gsl_monte_miser_free;
*gsl_monte_miser_params_get = *Math::GSL::Montec::gsl_monte_miser_params_get;
*gsl_monte_miser_params_set = *Math::GSL::Montec::gsl_monte_miser_params_set;
*gsl_monte_plain_integrate = *Math::GSL::Montec::gsl_monte_plain_integrate;
*gsl_monte_plain_alloc = *Math::GSL::Montec::gsl_monte_plain_alloc;
*gsl_monte_plain_init = *Math::GSL::Montec::gsl_monte_plain_init;
*gsl_monte_plain_free = *Math::GSL::Montec::gsl_monte_plain_free;
*gsl_monte_vegas_integrate = *Math::GSL::Montec::gsl_monte_vegas_integrate;
*gsl_monte_vegas_alloc = *Math::GSL::Montec::gsl_monte_vegas_alloc;
*gsl_monte_vegas_init = *Math::GSL::Montec::gsl_monte_vegas_init;
*gsl_monte_vegas_free = *Math::GSL::Montec::gsl_monte_vegas_free;
*gsl_monte_vegas_chisq = *Math::GSL::Montec::gsl_monte_vegas_chisq;
*gsl_monte_vegas_runval = *Math::GSL::Montec::gsl_monte_vegas_runval;
*gsl_monte_vegas_params_get = *Math::GSL::Montec::gsl_monte_vegas_params_get;
*gsl_monte_vegas_params_set = *Math::GSL::Montec::gsl_monte_vegas_params_set;
*gsl_error = *Math::GSL::Montec::gsl_error;
*gsl_stream_printf = *Math::GSL::Montec::gsl_stream_printf;
*gsl_strerror = *Math::GSL::Montec::gsl_strerror;
*gsl_set_error_handler = *Math::GSL::Montec::gsl_set_error_handler;
*gsl_set_error_handler_off = *Math::GSL::Montec::gsl_set_error_handler_off;
*gsl_set_stream_handler = *Math::GSL::Montec::gsl_set_stream_handler;
*gsl_set_stream = *Math::GSL::Montec::gsl_set_stream;

############# Class : Math::GSL::Monte::gsl_monte_function_struct ##############

package Math::GSL::Monte::gsl_monte_function_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Monte );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Montec::gsl_monte_function_struct_f_get;
*swig_f_set = *Math::GSL::Montec::gsl_monte_function_struct_f_set;
*swig_dim_get = *Math::GSL::Montec::gsl_monte_function_struct_dim_get;
*swig_dim_set = *Math::GSL::Montec::gsl_monte_function_struct_dim_set;
*swig_params_get = *Math::GSL::Montec::gsl_monte_function_struct_params_get;
*swig_params_set = *Math::GSL::Montec::gsl_monte_function_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Montec::new_gsl_monte_function_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Montec::delete_gsl_monte_function_struct($self);
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


############# Class : Math::GSL::Monte::gsl_monte_miser_state ##############

package Math::GSL::Monte::gsl_monte_miser_state;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Monte );
%OWNER = ();
%ITERATORS = ();
*swig_min_calls_get = *Math::GSL::Montec::gsl_monte_miser_state_min_calls_get;
*swig_min_calls_set = *Math::GSL::Montec::gsl_monte_miser_state_min_calls_set;
*swig_min_calls_per_bisection_get = *Math::GSL::Montec::gsl_monte_miser_state_min_calls_per_bisection_get;
*swig_min_calls_per_bisection_set = *Math::GSL::Montec::gsl_monte_miser_state_min_calls_per_bisection_set;
*swig_dither_get = *Math::GSL::Montec::gsl_monte_miser_state_dither_get;
*swig_dither_set = *Math::GSL::Montec::gsl_monte_miser_state_dither_set;
*swig_estimate_frac_get = *Math::GSL::Montec::gsl_monte_miser_state_estimate_frac_get;
*swig_estimate_frac_set = *Math::GSL::Montec::gsl_monte_miser_state_estimate_frac_set;
*swig_alpha_get = *Math::GSL::Montec::gsl_monte_miser_state_alpha_get;
*swig_alpha_set = *Math::GSL::Montec::gsl_monte_miser_state_alpha_set;
*swig_dim_get = *Math::GSL::Montec::gsl_monte_miser_state_dim_get;
*swig_dim_set = *Math::GSL::Montec::gsl_monte_miser_state_dim_set;
*swig_estimate_style_get = *Math::GSL::Montec::gsl_monte_miser_state_estimate_style_get;
*swig_estimate_style_set = *Math::GSL::Montec::gsl_monte_miser_state_estimate_style_set;
*swig_depth_get = *Math::GSL::Montec::gsl_monte_miser_state_depth_get;
*swig_depth_set = *Math::GSL::Montec::gsl_monte_miser_state_depth_set;
*swig_verbose_get = *Math::GSL::Montec::gsl_monte_miser_state_verbose_get;
*swig_verbose_set = *Math::GSL::Montec::gsl_monte_miser_state_verbose_set;
*swig_x_get = *Math::GSL::Montec::gsl_monte_miser_state_x_get;
*swig_x_set = *Math::GSL::Montec::gsl_monte_miser_state_x_set;
*swig_xmid_get = *Math::GSL::Montec::gsl_monte_miser_state_xmid_get;
*swig_xmid_set = *Math::GSL::Montec::gsl_monte_miser_state_xmid_set;
*swig_sigma_l_get = *Math::GSL::Montec::gsl_monte_miser_state_sigma_l_get;
*swig_sigma_l_set = *Math::GSL::Montec::gsl_monte_miser_state_sigma_l_set;
*swig_sigma_r_get = *Math::GSL::Montec::gsl_monte_miser_state_sigma_r_get;
*swig_sigma_r_set = *Math::GSL::Montec::gsl_monte_miser_state_sigma_r_set;
*swig_fmax_l_get = *Math::GSL::Montec::gsl_monte_miser_state_fmax_l_get;
*swig_fmax_l_set = *Math::GSL::Montec::gsl_monte_miser_state_fmax_l_set;
*swig_fmax_r_get = *Math::GSL::Montec::gsl_monte_miser_state_fmax_r_get;
*swig_fmax_r_set = *Math::GSL::Montec::gsl_monte_miser_state_fmax_r_set;
*swig_fmin_l_get = *Math::GSL::Montec::gsl_monte_miser_state_fmin_l_get;
*swig_fmin_l_set = *Math::GSL::Montec::gsl_monte_miser_state_fmin_l_set;
*swig_fmin_r_get = *Math::GSL::Montec::gsl_monte_miser_state_fmin_r_get;
*swig_fmin_r_set = *Math::GSL::Montec::gsl_monte_miser_state_fmin_r_set;
*swig_fsum_l_get = *Math::GSL::Montec::gsl_monte_miser_state_fsum_l_get;
*swig_fsum_l_set = *Math::GSL::Montec::gsl_monte_miser_state_fsum_l_set;
*swig_fsum_r_get = *Math::GSL::Montec::gsl_monte_miser_state_fsum_r_get;
*swig_fsum_r_set = *Math::GSL::Montec::gsl_monte_miser_state_fsum_r_set;
*swig_fsum2_l_get = *Math::GSL::Montec::gsl_monte_miser_state_fsum2_l_get;
*swig_fsum2_l_set = *Math::GSL::Montec::gsl_monte_miser_state_fsum2_l_set;
*swig_fsum2_r_get = *Math::GSL::Montec::gsl_monte_miser_state_fsum2_r_get;
*swig_fsum2_r_set = *Math::GSL::Montec::gsl_monte_miser_state_fsum2_r_set;
*swig_hits_l_get = *Math::GSL::Montec::gsl_monte_miser_state_hits_l_get;
*swig_hits_l_set = *Math::GSL::Montec::gsl_monte_miser_state_hits_l_set;
*swig_hits_r_get = *Math::GSL::Montec::gsl_monte_miser_state_hits_r_get;
*swig_hits_r_set = *Math::GSL::Montec::gsl_monte_miser_state_hits_r_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Montec::new_gsl_monte_miser_state(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Montec::delete_gsl_monte_miser_state($self);
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


############# Class : Math::GSL::Monte::gsl_monte_miser_params ##############

package Math::GSL::Monte::gsl_monte_miser_params;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Monte );
%OWNER = ();
%ITERATORS = ();
*swig_estimate_frac_get = *Math::GSL::Montec::gsl_monte_miser_params_estimate_frac_get;
*swig_estimate_frac_set = *Math::GSL::Montec::gsl_monte_miser_params_estimate_frac_set;
*swig_min_calls_get = *Math::GSL::Montec::gsl_monte_miser_params_min_calls_get;
*swig_min_calls_set = *Math::GSL::Montec::gsl_monte_miser_params_min_calls_set;
*swig_min_calls_per_bisection_get = *Math::GSL::Montec::gsl_monte_miser_params_min_calls_per_bisection_get;
*swig_min_calls_per_bisection_set = *Math::GSL::Montec::gsl_monte_miser_params_min_calls_per_bisection_set;
*swig_alpha_get = *Math::GSL::Montec::gsl_monte_miser_params_alpha_get;
*swig_alpha_set = *Math::GSL::Montec::gsl_monte_miser_params_alpha_set;
*swig_dither_get = *Math::GSL::Montec::gsl_monte_miser_params_dither_get;
*swig_dither_set = *Math::GSL::Montec::gsl_monte_miser_params_dither_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Montec::new_gsl_monte_miser_params(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Montec::delete_gsl_monte_miser_params($self);
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


############# Class : Math::GSL::Monte::gsl_monte_plain_state ##############

package Math::GSL::Monte::gsl_monte_plain_state;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Monte );
%OWNER = ();
%ITERATORS = ();
*swig_dim_get = *Math::GSL::Montec::gsl_monte_plain_state_dim_get;
*swig_dim_set = *Math::GSL::Montec::gsl_monte_plain_state_dim_set;
*swig_x_get = *Math::GSL::Montec::gsl_monte_plain_state_x_get;
*swig_x_set = *Math::GSL::Montec::gsl_monte_plain_state_x_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Montec::new_gsl_monte_plain_state(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Montec::delete_gsl_monte_plain_state($self);
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


############# Class : Math::GSL::Monte::gsl_monte_vegas_state ##############

package Math::GSL::Monte::gsl_monte_vegas_state;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Monte );
%OWNER = ();
%ITERATORS = ();
*swig_dim_get = *Math::GSL::Montec::gsl_monte_vegas_state_dim_get;
*swig_dim_set = *Math::GSL::Montec::gsl_monte_vegas_state_dim_set;
*swig_bins_max_get = *Math::GSL::Montec::gsl_monte_vegas_state_bins_max_get;
*swig_bins_max_set = *Math::GSL::Montec::gsl_monte_vegas_state_bins_max_set;
*swig_bins_get = *Math::GSL::Montec::gsl_monte_vegas_state_bins_get;
*swig_bins_set = *Math::GSL::Montec::gsl_monte_vegas_state_bins_set;
*swig_boxes_get = *Math::GSL::Montec::gsl_monte_vegas_state_boxes_get;
*swig_boxes_set = *Math::GSL::Montec::gsl_monte_vegas_state_boxes_set;
*swig_xi_get = *Math::GSL::Montec::gsl_monte_vegas_state_xi_get;
*swig_xi_set = *Math::GSL::Montec::gsl_monte_vegas_state_xi_set;
*swig_xin_get = *Math::GSL::Montec::gsl_monte_vegas_state_xin_get;
*swig_xin_set = *Math::GSL::Montec::gsl_monte_vegas_state_xin_set;
*swig_delx_get = *Math::GSL::Montec::gsl_monte_vegas_state_delx_get;
*swig_delx_set = *Math::GSL::Montec::gsl_monte_vegas_state_delx_set;
*swig_weight_get = *Math::GSL::Montec::gsl_monte_vegas_state_weight_get;
*swig_weight_set = *Math::GSL::Montec::gsl_monte_vegas_state_weight_set;
*swig_vol_get = *Math::GSL::Montec::gsl_monte_vegas_state_vol_get;
*swig_vol_set = *Math::GSL::Montec::gsl_monte_vegas_state_vol_set;
*swig_x_get = *Math::GSL::Montec::gsl_monte_vegas_state_x_get;
*swig_x_set = *Math::GSL::Montec::gsl_monte_vegas_state_x_set;
*swig_bin_get = *Math::GSL::Montec::gsl_monte_vegas_state_bin_get;
*swig_bin_set = *Math::GSL::Montec::gsl_monte_vegas_state_bin_set;
*swig_box_get = *Math::GSL::Montec::gsl_monte_vegas_state_box_get;
*swig_box_set = *Math::GSL::Montec::gsl_monte_vegas_state_box_set;
*swig_d_get = *Math::GSL::Montec::gsl_monte_vegas_state_d_get;
*swig_d_set = *Math::GSL::Montec::gsl_monte_vegas_state_d_set;
*swig_alpha_get = *Math::GSL::Montec::gsl_monte_vegas_state_alpha_get;
*swig_alpha_set = *Math::GSL::Montec::gsl_monte_vegas_state_alpha_set;
*swig_mode_get = *Math::GSL::Montec::gsl_monte_vegas_state_mode_get;
*swig_mode_set = *Math::GSL::Montec::gsl_monte_vegas_state_mode_set;
*swig_verbose_get = *Math::GSL::Montec::gsl_monte_vegas_state_verbose_get;
*swig_verbose_set = *Math::GSL::Montec::gsl_monte_vegas_state_verbose_set;
*swig_iterations_get = *Math::GSL::Montec::gsl_monte_vegas_state_iterations_get;
*swig_iterations_set = *Math::GSL::Montec::gsl_monte_vegas_state_iterations_set;
*swig_stage_get = *Math::GSL::Montec::gsl_monte_vegas_state_stage_get;
*swig_stage_set = *Math::GSL::Montec::gsl_monte_vegas_state_stage_set;
*swig_jac_get = *Math::GSL::Montec::gsl_monte_vegas_state_jac_get;
*swig_jac_set = *Math::GSL::Montec::gsl_monte_vegas_state_jac_set;
*swig_wtd_int_sum_get = *Math::GSL::Montec::gsl_monte_vegas_state_wtd_int_sum_get;
*swig_wtd_int_sum_set = *Math::GSL::Montec::gsl_monte_vegas_state_wtd_int_sum_set;
*swig_sum_wgts_get = *Math::GSL::Montec::gsl_monte_vegas_state_sum_wgts_get;
*swig_sum_wgts_set = *Math::GSL::Montec::gsl_monte_vegas_state_sum_wgts_set;
*swig_chi_sum_get = *Math::GSL::Montec::gsl_monte_vegas_state_chi_sum_get;
*swig_chi_sum_set = *Math::GSL::Montec::gsl_monte_vegas_state_chi_sum_set;
*swig_chisq_get = *Math::GSL::Montec::gsl_monte_vegas_state_chisq_get;
*swig_chisq_set = *Math::GSL::Montec::gsl_monte_vegas_state_chisq_set;
*swig_result_get = *Math::GSL::Montec::gsl_monte_vegas_state_result_get;
*swig_result_set = *Math::GSL::Montec::gsl_monte_vegas_state_result_set;
*swig_sigma_get = *Math::GSL::Montec::gsl_monte_vegas_state_sigma_get;
*swig_sigma_set = *Math::GSL::Montec::gsl_monte_vegas_state_sigma_set;
*swig_it_start_get = *Math::GSL::Montec::gsl_monte_vegas_state_it_start_get;
*swig_it_start_set = *Math::GSL::Montec::gsl_monte_vegas_state_it_start_set;
*swig_it_num_get = *Math::GSL::Montec::gsl_monte_vegas_state_it_num_get;
*swig_it_num_set = *Math::GSL::Montec::gsl_monte_vegas_state_it_num_set;
*swig_samples_get = *Math::GSL::Montec::gsl_monte_vegas_state_samples_get;
*swig_samples_set = *Math::GSL::Montec::gsl_monte_vegas_state_samples_set;
*swig_calls_per_box_get = *Math::GSL::Montec::gsl_monte_vegas_state_calls_per_box_get;
*swig_calls_per_box_set = *Math::GSL::Montec::gsl_monte_vegas_state_calls_per_box_set;
*swig_ostream_get = *Math::GSL::Montec::gsl_monte_vegas_state_ostream_get;
*swig_ostream_set = *Math::GSL::Montec::gsl_monte_vegas_state_ostream_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Montec::new_gsl_monte_vegas_state(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Montec::delete_gsl_monte_vegas_state($self);
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


############# Class : Math::GSL::Monte::gsl_monte_vegas_params ##############

package Math::GSL::Monte::gsl_monte_vegas_params;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Monte );
%OWNER = ();
%ITERATORS = ();
*swig_alpha_get = *Math::GSL::Montec::gsl_monte_vegas_params_alpha_get;
*swig_alpha_set = *Math::GSL::Montec::gsl_monte_vegas_params_alpha_set;
*swig_iterations_get = *Math::GSL::Montec::gsl_monte_vegas_params_iterations_get;
*swig_iterations_set = *Math::GSL::Montec::gsl_monte_vegas_params_iterations_set;
*swig_stage_get = *Math::GSL::Montec::gsl_monte_vegas_params_stage_get;
*swig_stage_set = *Math::GSL::Montec::gsl_monte_vegas_params_stage_set;
*swig_mode_get = *Math::GSL::Montec::gsl_monte_vegas_params_mode_get;
*swig_mode_set = *Math::GSL::Montec::gsl_monte_vegas_params_mode_set;
*swig_verbose_get = *Math::GSL::Montec::gsl_monte_vegas_params_verbose_get;
*swig_verbose_set = *Math::GSL::Montec::gsl_monte_vegas_params_verbose_set;
*swig_ostream_get = *Math::GSL::Montec::gsl_monte_vegas_params_ostream_get;
*swig_ostream_set = *Math::GSL::Montec::gsl_monte_vegas_params_ostream_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Montec::new_gsl_monte_vegas_params(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Montec::delete_gsl_monte_vegas_params($self);
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

package Math::GSL::Monte;

*GSL_MAJOR_VERSION = *Math::GSL::Montec::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Montec::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Montec::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Montec::GSL_NEGZERO;
*GSL_VEGAS_MODE_IMPORTANCE = *Math::GSL::Montec::GSL_VEGAS_MODE_IMPORTANCE;
*GSL_VEGAS_MODE_IMPORTANCE_ONLY = *Math::GSL::Montec::GSL_VEGAS_MODE_IMPORTANCE_ONLY;
*GSL_VEGAS_MODE_STRATIFIED = *Math::GSL::Montec::GSL_VEGAS_MODE_STRATIFIED;
*GSL_SUCCESS = *Math::GSL::Montec::GSL_SUCCESS;
*GSL_FAILURE = *Math::GSL::Montec::GSL_FAILURE;
*GSL_CONTINUE = *Math::GSL::Montec::GSL_CONTINUE;
*GSL_EDOM = *Math::GSL::Montec::GSL_EDOM;
*GSL_ERANGE = *Math::GSL::Montec::GSL_ERANGE;
*GSL_EFAULT = *Math::GSL::Montec::GSL_EFAULT;
*GSL_EINVAL = *Math::GSL::Montec::GSL_EINVAL;
*GSL_EFAILED = *Math::GSL::Montec::GSL_EFAILED;
*GSL_EFACTOR = *Math::GSL::Montec::GSL_EFACTOR;
*GSL_ESANITY = *Math::GSL::Montec::GSL_ESANITY;
*GSL_ENOMEM = *Math::GSL::Montec::GSL_ENOMEM;
*GSL_EBADFUNC = *Math::GSL::Montec::GSL_EBADFUNC;
*GSL_ERUNAWAY = *Math::GSL::Montec::GSL_ERUNAWAY;
*GSL_EMAXITER = *Math::GSL::Montec::GSL_EMAXITER;
*GSL_EZERODIV = *Math::GSL::Montec::GSL_EZERODIV;
*GSL_EBADTOL = *Math::GSL::Montec::GSL_EBADTOL;
*GSL_ETOL = *Math::GSL::Montec::GSL_ETOL;
*GSL_EUNDRFLW = *Math::GSL::Montec::GSL_EUNDRFLW;
*GSL_EOVRFLW = *Math::GSL::Montec::GSL_EOVRFLW;
*GSL_ELOSS = *Math::GSL::Montec::GSL_ELOSS;
*GSL_EROUND = *Math::GSL::Montec::GSL_EROUND;
*GSL_EBADLEN = *Math::GSL::Montec::GSL_EBADLEN;
*GSL_ENOTSQR = *Math::GSL::Montec::GSL_ENOTSQR;
*GSL_ESING = *Math::GSL::Montec::GSL_ESING;
*GSL_EDIVERGE = *Math::GSL::Montec::GSL_EDIVERGE;
*GSL_EUNSUP = *Math::GSL::Montec::GSL_EUNSUP;
*GSL_EUNIMPL = *Math::GSL::Montec::GSL_EUNIMPL;
*GSL_ECACHE = *Math::GSL::Montec::GSL_ECACHE;
*GSL_ETABLE = *Math::GSL::Montec::GSL_ETABLE;
*GSL_ENOPROG = *Math::GSL::Montec::GSL_ENOPROG;
*GSL_ENOPROGJ = *Math::GSL::Montec::GSL_ENOPROGJ;
*GSL_ETOLF = *Math::GSL::Montec::GSL_ETOLF;
*GSL_ETOLX = *Math::GSL::Montec::GSL_ETOLX;
*GSL_ETOLG = *Math::GSL::Montec::GSL_ETOLG;
*GSL_EOF = *Math::GSL::Montec::GSL_EOF;

@EXPORT_OK = qw/
               gsl_monte_miser_integrate 
               gsl_monte_miser_alloc 
               gsl_monte_miser_init 
               gsl_monte_miser_free 
               gsl_monte_plain_integrate 
               gsl_monte_plain_alloc 
               gsl_monte_plain_init 
               gsl_monte_plain_free 
               gsl_monte_vegas_integrate 
               gsl_monte_vegas_alloc 
               gsl_monte_vegas_init 
               gsl_monte_vegas_free 
               $GSL_VEGAS_MODE_IMPORTANCE 
               $GSL_VEGAS_MODE_IMPORTANCE_ONLY 
               $GSL_VEGAS_MODE_STRATIFIED 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Monte - Routines for multidimensional Monte Carlo integration

=head1 SYNOPSIS

This module is not yet implemented. Patches Welcome!

    use Math::GSL::Monte qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over 

=item * C<gsl_monte_miser_integrate >

=item * C<gsl_monte_miser_alloc >

=item * C<gsl_monte_miser_init >

=item * C<gsl_monte_miser_free >

=item * C<gsl_monte_plain_integrate >

=item * C<gsl_monte_plain_alloc >

=item * C<gsl_monte_plain_init >

=item * C<gsl_monte_plain_free >

=item * C<gsl_monte_vegas_integrate >

=item * C<gsl_monte_vegas_alloc >

=item * C<gsl_monte_vegas_init >

=item * C<gsl_monte_vegas_free >

=back

This module also includes the following constants :

=over

=item * $GSL_VEGAS_MODE_IMPORTANCE 

=item * $GSL_VEGAS_MODE_IMPORTANCE_ONLY 

=item * $GSL_VEGAS_MODE_STRATIFIED 

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want


=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut


1;

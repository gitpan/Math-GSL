# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Multifit;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Multifitc;
bootstrap Math::GSL::Multifit;
package Math::GSL::Multifit;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Multifit;

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

package Math::GSL::Multifit;

*gsl_multifit_linear_alloc = *Math::GSL::Multifitc::gsl_multifit_linear_alloc;
*gsl_multifit_linear_free = *Math::GSL::Multifitc::gsl_multifit_linear_free;
*gsl_multifit_linear = *Math::GSL::Multifitc::gsl_multifit_linear;
*gsl_multifit_linear_svd = *Math::GSL::Multifitc::gsl_multifit_linear_svd;
*gsl_multifit_wlinear = *Math::GSL::Multifitc::gsl_multifit_wlinear;
*gsl_multifit_wlinear_svd = *Math::GSL::Multifitc::gsl_multifit_wlinear_svd;
*gsl_multifit_linear_est = *Math::GSL::Multifitc::gsl_multifit_linear_est;
*gsl_multifit_linear_residuals = *Math::GSL::Multifitc::gsl_multifit_linear_residuals;
*gsl_multifit_gradient = *Math::GSL::Multifitc::gsl_multifit_gradient;
*gsl_multifit_covar = *Math::GSL::Multifitc::gsl_multifit_covar;
*gsl_multifit_fsolver_alloc = *Math::GSL::Multifitc::gsl_multifit_fsolver_alloc;
*gsl_multifit_fsolver_free = *Math::GSL::Multifitc::gsl_multifit_fsolver_free;
*gsl_multifit_fsolver_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_set;
*gsl_multifit_fsolver_iterate = *Math::GSL::Multifitc::gsl_multifit_fsolver_iterate;
*gsl_multifit_fsolver_name = *Math::GSL::Multifitc::gsl_multifit_fsolver_name;
*gsl_multifit_fsolver_position = *Math::GSL::Multifitc::gsl_multifit_fsolver_position;
*gsl_multifit_fdfsolver_alloc = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_alloc;
*gsl_multifit_fdfsolver_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_set;
*gsl_multifit_fdfsolver_iterate = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_iterate;
*gsl_multifit_fdfsolver_free = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_free;
*gsl_multifit_fdfsolver_name = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_name;
*gsl_multifit_fdfsolver_position = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_position;
*gsl_multifit_test_delta = *Math::GSL::Multifitc::gsl_multifit_test_delta;
*gsl_multifit_test_gradient = *Math::GSL::Multifitc::gsl_multifit_test_gradient;

############# Class : Math::GSL::Multifit::gsl_multifit_linear_workspace ##############

package Math::GSL::Multifit::gsl_multifit_linear_workspace;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multifit );
%OWNER = ();
%ITERATORS = ();
*swig_n_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_n_get;
*swig_n_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_n_set;
*swig_p_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_p_get;
*swig_p_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_p_set;
*swig_A_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_A_get;
*swig_A_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_A_set;
*swig_Q_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_Q_get;
*swig_Q_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_Q_set;
*swig_QSI_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_QSI_get;
*swig_QSI_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_QSI_set;
*swig_S_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_S_get;
*swig_S_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_S_set;
*swig_t_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_t_get;
*swig_t_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_t_set;
*swig_xt_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_xt_get;
*swig_xt_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_xt_set;
*swig_D_get = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_D_get;
*swig_D_set = *Math::GSL::Multifitc::gsl_multifit_linear_workspace_D_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multifitc::new_gsl_multifit_linear_workspace(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multifitc::delete_gsl_multifit_linear_workspace($self);
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


############# Class : Math::GSL::Multifit::gsl_multifit_function_struct ##############

package Math::GSL::Multifit::gsl_multifit_function_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multifit );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Multifitc::gsl_multifit_function_struct_f_get;
*swig_f_set = *Math::GSL::Multifitc::gsl_multifit_function_struct_f_set;
*swig_n_get = *Math::GSL::Multifitc::gsl_multifit_function_struct_n_get;
*swig_n_set = *Math::GSL::Multifitc::gsl_multifit_function_struct_n_set;
*swig_p_get = *Math::GSL::Multifitc::gsl_multifit_function_struct_p_get;
*swig_p_set = *Math::GSL::Multifitc::gsl_multifit_function_struct_p_set;
*swig_params_get = *Math::GSL::Multifitc::gsl_multifit_function_struct_params_get;
*swig_params_set = *Math::GSL::Multifitc::gsl_multifit_function_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multifitc::new_gsl_multifit_function_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multifitc::delete_gsl_multifit_function_struct($self);
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


############# Class : Math::GSL::Multifit::gsl_multifit_fsolver_type ##############

package Math::GSL::Multifit::gsl_multifit_fsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multifit );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_name_get;
*swig_name_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_name_set;
*swig_size_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_size_get;
*swig_size_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_size_set;
*swig_alloc_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_alloc_get;
*swig_alloc_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_alloc_set;
*swig_set_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_set_get;
*swig_set_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_iterate_set;
*swig_free_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_free_get;
*swig_free_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_free_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multifitc::new_gsl_multifit_fsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multifitc::delete_gsl_multifit_fsolver_type($self);
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


############# Class : Math::GSL::Multifit::gsl_multifit_fsolver ##############

package Math::GSL::Multifit::gsl_multifit_fsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multifit );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_get;
*swig_type_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_type_set;
*swig_function_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_function_get;
*swig_function_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_function_set;
*swig_x_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_x_get;
*swig_x_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_x_set;
*swig_f_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_f_get;
*swig_f_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_f_set;
*swig_dx_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_dx_get;
*swig_dx_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_dx_set;
*swig_state_get = *Math::GSL::Multifitc::gsl_multifit_fsolver_state_get;
*swig_state_set = *Math::GSL::Multifitc::gsl_multifit_fsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multifitc::new_gsl_multifit_fsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multifitc::delete_gsl_multifit_fsolver($self);
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


############# Class : Math::GSL::Multifit::gsl_multifit_function_fdf_struct ##############

package Math::GSL::Multifit::gsl_multifit_function_fdf_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multifit );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_f_get;
*swig_f_set = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_f_set;
*swig_df_get = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_df_get;
*swig_df_set = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_df_set;
*swig_fdf_get = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_fdf_get;
*swig_fdf_set = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_fdf_set;
*swig_n_get = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_n_get;
*swig_n_set = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_n_set;
*swig_p_get = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_p_get;
*swig_p_set = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_p_set;
*swig_params_get = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_params_get;
*swig_params_set = *Math::GSL::Multifitc::gsl_multifit_function_fdf_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multifitc::new_gsl_multifit_function_fdf_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multifitc::delete_gsl_multifit_function_fdf_struct($self);
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


############# Class : Math::GSL::Multifit::gsl_multifit_fdfsolver_type ##############

package Math::GSL::Multifit::gsl_multifit_fdfsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multifit );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_name_get;
*swig_name_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_name_set;
*swig_size_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_size_get;
*swig_size_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_size_set;
*swig_alloc_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_alloc_get;
*swig_alloc_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_alloc_set;
*swig_set_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_set_get;
*swig_set_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_iterate_set;
*swig_free_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_free_get;
*swig_free_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_free_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multifitc::new_gsl_multifit_fdfsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multifitc::delete_gsl_multifit_fdfsolver_type($self);
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


############# Class : Math::GSL::Multifit::gsl_multifit_fdfsolver ##############

package Math::GSL::Multifit::gsl_multifit_fdfsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multifit );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_get;
*swig_type_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_type_set;
*swig_fdf_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_fdf_get;
*swig_fdf_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_fdf_set;
*swig_x_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_x_get;
*swig_x_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_x_set;
*swig_f_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_f_get;
*swig_f_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_f_set;
*swig_J_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_J_get;
*swig_J_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_J_set;
*swig_dx_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_dx_get;
*swig_dx_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_dx_set;
*swig_state_get = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_state_get;
*swig_state_set = *Math::GSL::Multifitc::gsl_multifit_fdfsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multifitc::new_gsl_multifit_fdfsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multifitc::delete_gsl_multifit_fdfsolver($self);
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

package Math::GSL::Multifit;

*GSL_MAJOR_VERSION = *Math::GSL::Multifitc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Multifitc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Multifitc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Multifitc::GSL_NEGZERO;

my %__gsl_multifit_fdfsolver_lmder_hash;
tie %__gsl_multifit_fdfsolver_lmder_hash,"Math::GSL::Multifit::gsl_multifit_fdfsolver_type", $Math::GSL::Multifitc::gsl_multifit_fdfsolver_lmder;
$gsl_multifit_fdfsolver_lmder= \%__gsl_multifit_fdfsolver_lmder_hash;
bless $gsl_multifit_fdfsolver_lmder, Math::GSL::Multifit::gsl_multifit_fdfsolver_type;

my %__gsl_multifit_fdfsolver_lmsder_hash;
tie %__gsl_multifit_fdfsolver_lmsder_hash,"Math::GSL::Multifit::gsl_multifit_fdfsolver_type", $Math::GSL::Multifitc::gsl_multifit_fdfsolver_lmsder;
$gsl_multifit_fdfsolver_lmsder= \%__gsl_multifit_fdfsolver_lmsder_hash;
bless $gsl_multifit_fdfsolver_lmsder, Math::GSL::Multifit::gsl_multifit_fdfsolver_type;

@EXPORT_OK = qw/
               gsl_multifit_linear_alloc 
               gsl_multifit_linear_free 
               gsl_multifit_linear 
               gsl_multifit_linear_svd 
               gsl_multifit_wlinear 
               gsl_multifit_wlinear_svd 
               gsl_multifit_linear_est 
               gsl_multifit_linear_residuals 
               gsl_multifit_gradient 
               gsl_multifit_covar 
               gsl_multifit_fsolver_alloc 
               gsl_multifit_fsolver_free 
               gsl_multifit_fsolver_set 
               gsl_multifit_fsolver_iterate 
               gsl_multifit_fsolver_name 
               gsl_multifit_fsolver_position 
               gsl_multifit_fdfsolver_alloc 
               gsl_multifit_fdfsolver_set 
               gsl_multifit_fdfsolver_iterate 
               gsl_multifit_fdfsolver_free 
               gsl_multifit_fdfsolver_name 
               gsl_multifit_fdfsolver_position 
               gsl_multifit_test_delta 
               gsl_multifit_test_gradient
               $gsl_multifit_fdfsolver_lmder
               $gsl_multifit_fdfsolver_lmsder; 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Multifit - Least-squares functions for a general linear model with multiple parameters

=head1 SYNOPSIS

use Math::GSL::Multifit qw /:all/;

=head1 DESCRIPTION

The functions in this module perform least-squares fits to a general linear model, y = X c where y is a vector of n observations, X is an n by p matrix of predictor variables, and the elements of the vector c are the p unknown best-fit parameters which are to be estimated.

Here is a list of all the functions in this module :

=over 

=item C<gsl_multifit_linear_alloc($n, $p)> - This function allocates a workspace for fitting a model to $n observations using $p parameters. 

=item C<gsl_multifit_linear_free($work)> - This function frees the memory associated with the workspace w. 

=item C<gsl_multifit_linear($X, $y, $c, $cov, $work)> - This function computes the best-fit parameters vector $c of the model y = X c for the observations vector $y and the matrix of predictor variables $X. The variance-covariance matrix of the model parameters vector $cov is estimated from the scatter of the observations about the best-fit. The sum of squares of the residuals from the best-fit, \chi^2, is returned after 0 if the operation succeeded, 1 otherwise. If the coefficient of determination is desired, it can be computed from the expression R^2 = 1 - \chi^2 / TSS, where the total sum of squares (TSS) of the observations y may be computed from gsl_stats_tss. The best-fit is found by singular value decomposition of the matrix $X using the preallocated workspace provided in $work. The modified Golub-Reinsch SVD algorithm is used, with column scaling to improve the accuracy of the singular values. Any components which have zero singular value (to machine precision) are discarded from the fit.

=item C<gsl_multifit_linear_svd($X, $y, $tol, $c, $cov, $work)> - This function computes the best-fit parameters c of the model y = X c for the observations vector $y and the matrix of predictor variables $X. The variance-covariance matrix of the model parameters vector $cov is estimated from the scatter of the observations about the best-fit. The sum of squares of the residuals from the best-fit, \chi^2, is returned after 0 if the operation succeeded, 1 otherwise. If the coefficient of determination is desired, it can be computed from the expression R^2 = 1 - \chi^2 / TSS, where the total sum of squares (TSS) of the observations y may be computed from gsl_stats_tss. In this second form of the function the components are discarded if the ratio of singular values s_i/s_0 falls below the user-specified tolerance $tol, and the effective rank is returned after the sum of squares of the residuals from the best-fit.

=item C<gsl_multifit_wlinear($X, $w, $y, $c, $cov, $work> - This function computes the best-fit parameters vector $c of the weighted model y = X c for the observations y with weights $w and the matrix of predictor variables $X. The covariance matrix of the model parameters $cov is computed with the given weights. The weighted sum of squares of the residuals from the best-fit, \chi^2, is returned after 0 if the operation succeeded, 1 otherwise. If the coefficient of determination is desired, it can be computed from the expression R^2 = 1 - \chi^2 / WTSS, where the weighted total sum of squares (WTSS) of the observations y may be computed from gsl_stats_wtss. The best-fit is found by singular value decomposition of the matrix $X using the preallocated workspace provided in $work. Any components which have zero singular value (to machine precision) are discarded from the fit.

=item C<gsl_multifit_wlinear_svd($X, $w, $y, $tol, $rank, $c, $cov, $work) > This function computes the best-fit parameters vector $c of the weighted model y = X c for the observations y with weights $w and the matrix of predictor variables $X. The covariance matrix of the model parameters $cov is computed with the given weights. The weighted sum of squares of the residuals from the best-fit, \chi^2, is returned after 0 if the operation succeeded, 1 otherwise. If the coefficient of determination is desired, it can be computed from the expression R^2 = 1 - \chi^2 / WTSS, where the weighted total sum of squares (WTSS) of the observations y may be computed from gsl_stats_wtss. The best-fit is found by singular value decomposition of the matrix $X using the preallocated workspace provided in $work. In this second form of the function the components are discarded if the ratio of singular values s_i/s_0 falls below the user-specified tolerance $tol, and the effective rank is returned after the sum of squares of the residuals from the best-fit.. 

=item C<gsl_multifit_linear_est($x, $c, $cov)> - This function uses the best-fit multilinear regression coefficients vector $c and their covariance matrix $cov to compute the fitted function value $y and its standard deviation $y_err for the model y = x.c at the point $x, in the form of a vector. The functions returns 3 values in this order : 0 if the operation succeeded, 1 otherwise, the fittes function value and its standard deviation. 

=item C<gsl_multifit_linear_residuals($X, $y, $c, $r)> - This function computes the vector of residuals r = y - X c for the observations vector $y, coefficients vector $c and matrix of predictor variables $X. $r is also a vector.

=item C<gsl_multifit_gradient($J, $f, $g)> - This function computes the gradient $g of \Phi(x) = (1/2) ||F(x)||^2 from the Jacobian matrix $J and the function values $f, using the formula $g = $J^T $f. $g and $f are vectors.

=item C<gsl_multifit_test_gradient($g, $epsabas)> - This function tests the residual gradient vector $g against the absolute error bound $epsabs. Mathematically, the gradient should be exactly zero at the minimum. The test returns $GSL_SUCCESS if the following condition is achieved, \sum_i |g_i| < $epsabs and returns $GSL_CONTINUE otherwise. This criterion is suitable for situations where the precise location of the minimum, x, is unimportant provided a value can be found where the gradient is small enough. 

=item C<gsl_multifit_test_delta($dx, $x, $epsabs, $epsrel)> - This function tests for the convergence of the sequence by comparing the last step vector $dx with the absolute error $epsabs and relative error $epsrel to the current position x. The test returns $GSL_SUCCESS if the following condition is achieved, |dx_i| < epsabs + epsrel |x_i| for each component of x and returns $GSL_CONTINUE otherwise. 

=back

The following functions are not yet implemented. Patches Welcome! 

=over

=item C<gsl_multifit_covar >

=item C<gsl_multifit_fsolver_alloc($T, $n, $p)>

=item C<gsl_multifit_fsolver_free >

=item C<gsl_multifit_fsolver_set >

=item C<gsl_multifit_fsolver_iterate >

=item C<gsl_multifit_fsolver_name >

=item C<gsl_multifit_fsolver_position >

=item C<gsl_multifit_fdfsolver_alloc >

=item C<gsl_multifit_fdfsolver_set >

=item C<gsl_multifit_fdfsolver_iterate >

=item C<gsl_multifit_fdfsolver_free >

=item C<gsl_multifit_fdfsolver_name >

=item C<gsl_multifit_fdfsolver_position >


=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

 Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want

=head1 EXAMPLES



=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::ODEIV;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::ODEIVc;
bootstrap Math::GSL::ODEIV;
package Math::GSL::ODEIV;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::ODEIV;

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

package Math::GSL::ODEIV;

*gsl_odeiv_step_alloc = *Math::GSL::ODEIVc::gsl_odeiv_step_alloc;
*gsl_odeiv_step_reset = *Math::GSL::ODEIVc::gsl_odeiv_step_reset;
*gsl_odeiv_step_free = *Math::GSL::ODEIVc::gsl_odeiv_step_free;
*gsl_odeiv_step_name = *Math::GSL::ODEIVc::gsl_odeiv_step_name;
*gsl_odeiv_step_order = *Math::GSL::ODEIVc::gsl_odeiv_step_order;
*gsl_odeiv_step_apply = *Math::GSL::ODEIVc::gsl_odeiv_step_apply;
*gsl_odeiv_control_alloc = *Math::GSL::ODEIVc::gsl_odeiv_control_alloc;
*gsl_odeiv_control_init = *Math::GSL::ODEIVc::gsl_odeiv_control_init;
*gsl_odeiv_control_free = *Math::GSL::ODEIVc::gsl_odeiv_control_free;
*gsl_odeiv_control_hadjust = *Math::GSL::ODEIVc::gsl_odeiv_control_hadjust;
*gsl_odeiv_control_name = *Math::GSL::ODEIVc::gsl_odeiv_control_name;
*gsl_odeiv_control_standard_new = *Math::GSL::ODEIVc::gsl_odeiv_control_standard_new;
*gsl_odeiv_control_y_new = *Math::GSL::ODEIVc::gsl_odeiv_control_y_new;
*gsl_odeiv_control_yp_new = *Math::GSL::ODEIVc::gsl_odeiv_control_yp_new;
*gsl_odeiv_control_scaled_new = *Math::GSL::ODEIVc::gsl_odeiv_control_scaled_new;
*gsl_odeiv_evolve_alloc = *Math::GSL::ODEIVc::gsl_odeiv_evolve_alloc;
*gsl_odeiv_evolve_apply = *Math::GSL::ODEIVc::gsl_odeiv_evolve_apply;
*gsl_odeiv_evolve_reset = *Math::GSL::ODEIVc::gsl_odeiv_evolve_reset;
*gsl_odeiv_evolve_free = *Math::GSL::ODEIVc::gsl_odeiv_evolve_free;

############# Class : Math::GSL::ODEIV::gsl_odeiv_system ##############

package Math::GSL::ODEIV::gsl_odeiv_system;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::ODEIV );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::ODEIVc::gsl_odeiv_system_function_get;
*swig_function_set = *Math::GSL::ODEIVc::gsl_odeiv_system_function_set;
*swig_jacobian_get = *Math::GSL::ODEIVc::gsl_odeiv_system_jacobian_get;
*swig_jacobian_set = *Math::GSL::ODEIVc::gsl_odeiv_system_jacobian_set;
*swig_dimension_get = *Math::GSL::ODEIVc::gsl_odeiv_system_dimension_get;
*swig_dimension_set = *Math::GSL::ODEIVc::gsl_odeiv_system_dimension_set;
*swig_params_get = *Math::GSL::ODEIVc::gsl_odeiv_system_params_get;
*swig_params_set = *Math::GSL::ODEIVc::gsl_odeiv_system_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::ODEIVc::new_gsl_odeiv_system(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::ODEIVc::delete_gsl_odeiv_system($self);
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


############# Class : Math::GSL::ODEIV::gsl_odeiv_step_type ##############

package Math::GSL::ODEIV::gsl_odeiv_step_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::ODEIV );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_name_get;
*swig_name_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_name_set;
*swig_can_use_dydt_in_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_can_use_dydt_in_get;
*swig_can_use_dydt_in_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_can_use_dydt_in_set;
*swig_gives_exact_dydt_out_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_gives_exact_dydt_out_get;
*swig_gives_exact_dydt_out_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_gives_exact_dydt_out_set;
*swig_alloc_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_alloc_get;
*swig_alloc_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_alloc_set;
*swig_apply_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_apply_get;
*swig_apply_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_apply_set;
*swig_reset_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_reset_get;
*swig_reset_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_reset_set;
*swig_order_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_order_get;
*swig_order_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_order_set;
*swig_free_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_free_get;
*swig_free_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_free_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::ODEIVc::new_gsl_odeiv_step_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::ODEIVc::delete_gsl_odeiv_step_type($self);
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


############# Class : Math::GSL::ODEIV::gsl_odeiv_step ##############

package Math::GSL::ODEIV::gsl_odeiv_step;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::ODEIV );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::ODEIVc::gsl_odeiv_step_type_get;
*swig_type_set = *Math::GSL::ODEIVc::gsl_odeiv_step_type_set;
*swig_dimension_get = *Math::GSL::ODEIVc::gsl_odeiv_step_dimension_get;
*swig_dimension_set = *Math::GSL::ODEIVc::gsl_odeiv_step_dimension_set;
*swig_state_get = *Math::GSL::ODEIVc::gsl_odeiv_step_state_get;
*swig_state_set = *Math::GSL::ODEIVc::gsl_odeiv_step_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::ODEIVc::new_gsl_odeiv_step(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::ODEIVc::delete_gsl_odeiv_step($self);
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


############# Class : Math::GSL::ODEIV::gsl_odeiv_control_type ##############

package Math::GSL::ODEIV::gsl_odeiv_control_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::ODEIV );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::ODEIVc::gsl_odeiv_control_type_name_get;
*swig_name_set = *Math::GSL::ODEIVc::gsl_odeiv_control_type_name_set;
*swig_alloc_get = *Math::GSL::ODEIVc::gsl_odeiv_control_type_alloc_get;
*swig_alloc_set = *Math::GSL::ODEIVc::gsl_odeiv_control_type_alloc_set;
*swig_init_get = *Math::GSL::ODEIVc::gsl_odeiv_control_type_init_get;
*swig_init_set = *Math::GSL::ODEIVc::gsl_odeiv_control_type_init_set;
*swig_hadjust_get = *Math::GSL::ODEIVc::gsl_odeiv_control_type_hadjust_get;
*swig_hadjust_set = *Math::GSL::ODEIVc::gsl_odeiv_control_type_hadjust_set;
*swig_free_get = *Math::GSL::ODEIVc::gsl_odeiv_control_type_free_get;
*swig_free_set = *Math::GSL::ODEIVc::gsl_odeiv_control_type_free_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::ODEIVc::new_gsl_odeiv_control_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::ODEIVc::delete_gsl_odeiv_control_type($self);
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


############# Class : Math::GSL::ODEIV::gsl_odeiv_control ##############

package Math::GSL::ODEIV::gsl_odeiv_control;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::ODEIV );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::ODEIVc::gsl_odeiv_control_type_get;
*swig_type_set = *Math::GSL::ODEIVc::gsl_odeiv_control_type_set;
*swig_state_get = *Math::GSL::ODEIVc::gsl_odeiv_control_state_get;
*swig_state_set = *Math::GSL::ODEIVc::gsl_odeiv_control_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::ODEIVc::new_gsl_odeiv_control(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::ODEIVc::delete_gsl_odeiv_control($self);
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


############# Class : Math::GSL::ODEIV::gsl_odeiv_evolve ##############

package Math::GSL::ODEIV::gsl_odeiv_evolve;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::ODEIV );
%OWNER = ();
%ITERATORS = ();
*swig_dimension_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_dimension_get;
*swig_dimension_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_dimension_set;
*swig_y0_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_y0_get;
*swig_y0_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_y0_set;
*swig_yerr_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_yerr_get;
*swig_yerr_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_yerr_set;
*swig_dydt_in_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_dydt_in_get;
*swig_dydt_in_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_dydt_in_set;
*swig_dydt_out_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_dydt_out_get;
*swig_dydt_out_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_dydt_out_set;
*swig_last_step_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_last_step_get;
*swig_last_step_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_last_step_set;
*swig_count_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_count_get;
*swig_count_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_count_set;
*swig_failed_steps_get = *Math::GSL::ODEIVc::gsl_odeiv_evolve_failed_steps_get;
*swig_failed_steps_set = *Math::GSL::ODEIVc::gsl_odeiv_evolve_failed_steps_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::ODEIVc::new_gsl_odeiv_evolve(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::ODEIVc::delete_gsl_odeiv_evolve($self);
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

package Math::GSL::ODEIV;

*GSL_MAJOR_VERSION = *Math::GSL::ODEIVc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::ODEIVc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::ODEIVc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::ODEIVc::GSL_NEGZERO;

my %__gsl_odeiv_step_rk2_hash;
tie %__gsl_odeiv_step_rk2_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rk2;
$gsl_odeiv_step_rk2= \%__gsl_odeiv_step_rk2_hash;
bless $gsl_odeiv_step_rk2, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_rk4_hash;
tie %__gsl_odeiv_step_rk4_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rk4;
$gsl_odeiv_step_rk4= \%__gsl_odeiv_step_rk4_hash;
bless $gsl_odeiv_step_rk4, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_rkf45_hash;
tie %__gsl_odeiv_step_rkf45_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rkf45;
$gsl_odeiv_step_rkf45= \%__gsl_odeiv_step_rkf45_hash;
bless $gsl_odeiv_step_rkf45, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_rkck_hash;
tie %__gsl_odeiv_step_rkck_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rkck;
$gsl_odeiv_step_rkck= \%__gsl_odeiv_step_rkck_hash;
bless $gsl_odeiv_step_rkck, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_rk8pd_hash;
tie %__gsl_odeiv_step_rk8pd_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rk8pd;
$gsl_odeiv_step_rk8pd= \%__gsl_odeiv_step_rk8pd_hash;
bless $gsl_odeiv_step_rk8pd, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_rk2imp_hash;
tie %__gsl_odeiv_step_rk2imp_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rk2imp;
$gsl_odeiv_step_rk2imp= \%__gsl_odeiv_step_rk2imp_hash;
bless $gsl_odeiv_step_rk2imp, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_rk2simp_hash;
tie %__gsl_odeiv_step_rk2simp_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rk2simp;
$gsl_odeiv_step_rk2simp= \%__gsl_odeiv_step_rk2simp_hash;
bless $gsl_odeiv_step_rk2simp, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_rk4imp_hash;
tie %__gsl_odeiv_step_rk4imp_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_rk4imp;
$gsl_odeiv_step_rk4imp= \%__gsl_odeiv_step_rk4imp_hash;
bless $gsl_odeiv_step_rk4imp, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_bsimp_hash;
tie %__gsl_odeiv_step_bsimp_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_bsimp;
$gsl_odeiv_step_bsimp= \%__gsl_odeiv_step_bsimp_hash;
bless $gsl_odeiv_step_bsimp, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_gear1_hash;
tie %__gsl_odeiv_step_gear1_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_gear1;
$gsl_odeiv_step_gear1= \%__gsl_odeiv_step_gear1_hash;
bless $gsl_odeiv_step_gear1, Math::GSL::ODEIV::gsl_odeiv_step_type;

my %__gsl_odeiv_step_gear2_hash;
tie %__gsl_odeiv_step_gear2_hash,"Math::GSL::ODEIV::gsl_odeiv_step_type", $Math::GSL::ODEIVc::gsl_odeiv_step_gear2;
$gsl_odeiv_step_gear2= \%__gsl_odeiv_step_gear2_hash;
bless $gsl_odeiv_step_gear2, Math::GSL::ODEIV::gsl_odeiv_step_type;
*GSL_ODEIV_HADJ_INC = *Math::GSL::ODEIVc::GSL_ODEIV_HADJ_INC;
*GSL_ODEIV_HADJ_NIL = *Math::GSL::ODEIVc::GSL_ODEIV_HADJ_NIL;
*GSL_ODEIV_HADJ_DEC = *Math::GSL::ODEIVc::GSL_ODEIV_HADJ_DEC;

@EXPORT_OK = qw/
               gsl_odeiv_step_alloc 
               gsl_odeiv_step_reset 
               gsl_odeiv_step_free 
               gsl_odeiv_step_name 
               gsl_odeiv_step_order 
               gsl_odeiv_step_apply 
               gsl_odeiv_control_alloc 
               gsl_odeiv_control_init 
               gsl_odeiv_control_free 
               gsl_odeiv_control_hadjust 
               gsl_odeiv_control_name 
               gsl_odeiv_control_standard_new 
               gsl_odeiv_control_y_new 
               gsl_odeiv_control_yp_new 
               gsl_odeiv_control_scaled_new 
               gsl_odeiv_evolve_alloc 
               gsl_odeiv_evolve_apply 
               gsl_odeiv_evolve_reset 
               gsl_odeiv_evolve_free 
               $gsl_odeiv_step_rk2
               $gsl_odeiv_step_rk4
               $gsl_odeiv_step_rkf45
               $gsl_odeiv_step_rkck
               $gsl_odeiv_step_rk8pd
               $gsl_odeiv_step_rk2imp
               $gsl_odeiv_step_rk2simp
               $gsl_odeiv_step_rk4imp
               $gsl_odeiv_step_bsimp
               $gsl_odeiv_step_gear1
               $gsl_odeiv_step_gear2
               $GSL_ODEIV_HADJ_INC 
               $GSL_ODEIV_HADJ_NIL 
               $GSL_ODEIV_HADJ_DEC
               $gsl_odeiv_control_standard             
	/;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::ODEIV - functions for solving ordinary differential equation (ODE) initial value problems

=head1 SYNOPSIS

use Math::GSL::ODEIV qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over

=item * C<gsl_odeiv_step_alloc($T, $dim)> - This function returns a pointer to a newly allocated instance of a stepping function of type $T for a system of $dim dimensions.$T must be one of the step type constant above. 

=item * C<gsl_odeiv_step_reset($s)> - This function resets the stepping function $s. It should be used whenever the next use of s will not be a continuation of a previous step.

=item * C<gsl_odeiv_step_free($s)> - This function frees all the memory associated with the stepping function $s.

=item * C<gsl_odeiv_step_name($s)> - This function returns a pointer to the name of the stepping function.

=item * C<gsl_odeiv_step_order($s)> - This function returns the order of the stepping function on the previous step. This order can vary if the stepping function itself is adaptive.

=item * C<gsl_odeiv_step_apply >

=item * C<gsl_odeiv_control_alloc($T)> - This function returns a pointer to a newly allocated instance of a control function of type $T. This function is only needed for defining new types of control functions. For most purposes the standard control functions described above should be sufficient. $T is a gsl_odeiv_control_type. 

=item * C<gsl_odeiv_control_init($c, $eps_abs, $eps_rel, $a_y, $a_dydt) > - This function initializes the control function c with the parameters eps_abs (absolute error), eps_rel (relative error), a_y (scaling factor for y) and a_dydt (scaling factor for derivatives). 

=item * C<gsl_odeiv_control_free >

=item * C<gsl_odeiv_control_hadjust >

=item * C<gsl_odeiv_control_name >

=item * C<gsl_odeiv_control_standard_new($eps_abs, $eps_rel, $a_y, $a_dydt)> - The standard control object is a four parameter heuristic based on absolute and relative errors $eps_abs and $eps_rel, and scaling factors $a_y and $a_dydt for the system state y(t) and derivatives y'(t) respectively. The step-size adjustment procedure for this method begins by computing the desired error level D_i for each component, D_i = eps_abs + eps_rel * (a_y |y_i| + a_dydt h |y'_i|) and comparing it with the observed error E_i = |yerr_i|. If the observed error E exceeds the desired error level D by more than 10% for any component then the method reduces the step-size by an appropriate factor, h_new = h_old * S * (E/D)^(-1/q) where q is the consistency order of the method (e.g. q=4 for 4(5) embedded RK), and S is a safety factor of 0.9. The ratio E/D is taken to be the maximum of the ratios E_i/D_i. If the observed error E is less than 50% of the desired error level D for the maximum ratio E_i/D_i then the algorithm takes the opportunity to increase the step-size to bring the error in line with the desired level, h_new = h_old * S * (E/D)^(-1/(q+1)) This encompasses all the standard error scaling methods. To avoid uncontrolled changes in the stepsize, the overall scaling factor is limited to the range 1/5 to 5. 

=item * C<gsl_odeiv_control_y_new($eps_abs, $eps_rel)> - This function creates a new control object which will keep the local error on each step within an absolute error of $eps_abs and relative error of $eps_rel with respect to the solution y_i(t). This is equivalent to the standard control object with a_y=1 and a_dydt=0. 

=item * C<gsl_odeiv_control_yp_new($eps_abs, $eps_rel)> - This function creates a new control object which will keep the local error on each step within an absolute error of $eps_abs and relative error of $eps_rel with respect to the derivatives of the solution y'_i(t). This is equivalent to the standard control object with a_y=0 and a_dydt=1. 

=item * C<gsl_odeiv_control_scaled_new($eps_abs, $eps_rel, $a_y, $a_dydt, $scale_abs, $dim) > - This function creates a new control object which uses the same algorithm as gsl_odeiv_control_standard_new but with an absolute error which is scaled for each component by the array reference $scale_abs. The formula for D_i for this control object is, D_i = eps_abs * s_i + eps_rel * (a_y |y_i| + a_dydt h |y'_i|) where s_i is the i-th component of the array scale_abs. The same error control heuristic is used by the Matlab ode suite. 

=item * C<gsl_odeiv_evolve_alloc($dim)> - This function returns a pointer to a newly allocated instance of an evolution function for a system of $dim dimensions.

=item * C<gsl_odeiv_evolve_apply >

=item * C<gsl_odeiv_evolve_reset($e)> - This function resets the evolution function $e. It should be used whenever the next use of $e will not be a continuation of a previous step.

=item * C<gsl_odeiv_evolve_free($e)> - This function frees all the memory associated with the evolution function $e. 

=back

This module also includes the following constants :

=over

=item * C<$GSL_ODEIV_HADJ_INC>

=item * C<$GSL_ODEIV_HADJ_NIL>

=item * C<$GSL_ODEIV_HADJ_DEC>

=back

=head2 Step Type

=over

=item * C<$gsl_odeiv_step_rk2> - Embedded Runge-Kutta (2, 3) method.

=item * C<$gsl_odeiv_step_rk4> - 4th order (classical) Runge-Kutta. The error estimate is obtained by halving the step-size. For more efficient estimate of the error, use the Runge-Kutta-Fehlberg method described below.

=item * C<$gsl_odeiv_step_rkf45> - Embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good general-purpose integrator. 

=item * C<$gsl_odeiv_step_rkck> - Embedded Runge-Kutta Cash-Karp (4, 5) method.

=item * C<$gsl_odeiv_step_rk8pd> - Embedded Runge-Kutta Prince-Dormand (8,9) method. 

=item * C<$gsl_odeiv_step_rk2imp> - Implicit 2nd order Runge-Kutta at Gaussian points. 

=item * C<$gsl_odeiv_step_rk2simp>

=item * C<$gsl_odeiv_step_rk4imp> - Implicit 4th order Runge-Kutta at Gaussian points. 

=item * C<$gsl_odeiv_step_bsimp> - Implicit Bulirsch-Stoer method of Bader and Deuflhard. This algorithm requires the Jacobian. 

=item * C<$gsl_odeiv_step_gear1> - M=1 implicit Gear method. 

=item * C<$gsl_odeiv_step_gear2> - M=2 implicit Gear method.

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

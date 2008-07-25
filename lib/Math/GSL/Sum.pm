# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.36
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Sum;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Sumc;
bootstrap Math::GSL::Sum;
package Math::GSL::Sum;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Sum;

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

package Math::GSL::Sum;

*gsl_sum_levin_u_alloc = *Math::GSL::Sumc::gsl_sum_levin_u_alloc;
*gsl_sum_levin_u_free = *Math::GSL::Sumc::gsl_sum_levin_u_free;
*gsl_sum_levin_u_accel = *Math::GSL::Sumc::gsl_sum_levin_u_accel;
*gsl_sum_levin_u_minmax = *Math::GSL::Sumc::gsl_sum_levin_u_minmax;
*gsl_sum_levin_u_step = *Math::GSL::Sumc::gsl_sum_levin_u_step;
*gsl_sum_levin_utrunc_alloc = *Math::GSL::Sumc::gsl_sum_levin_utrunc_alloc;
*gsl_sum_levin_utrunc_free = *Math::GSL::Sumc::gsl_sum_levin_utrunc_free;
*gsl_sum_levin_utrunc_accel = *Math::GSL::Sumc::gsl_sum_levin_utrunc_accel;
*gsl_sum_levin_utrunc_minmax = *Math::GSL::Sumc::gsl_sum_levin_utrunc_minmax;
*gsl_sum_levin_utrunc_step = *Math::GSL::Sumc::gsl_sum_levin_utrunc_step;

############# Class : Math::GSL::Sum::gsl_sum_levin_u_workspace ##############

package Math::GSL::Sum::gsl_sum_levin_u_workspace;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Sum );
%OWNER = ();
%ITERATORS = ();
*swig_size_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_size_get;
*swig_size_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_size_set;
*swig_i_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_i_get;
*swig_i_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_i_set;
*swig_terms_used_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_terms_used_get;
*swig_terms_used_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_terms_used_set;
*swig_sum_plain_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_sum_plain_get;
*swig_sum_plain_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_sum_plain_set;
*swig_q_num_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_q_num_get;
*swig_q_num_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_q_num_set;
*swig_q_den_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_q_den_get;
*swig_q_den_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_q_den_set;
*swig_dq_num_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_dq_num_get;
*swig_dq_num_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_dq_num_set;
*swig_dq_den_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_dq_den_get;
*swig_dq_den_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_dq_den_set;
*swig_dsum_get = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_dsum_get;
*swig_dsum_set = *Math::GSL::Sumc::gsl_sum_levin_u_workspace_dsum_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Sumc::new_gsl_sum_levin_u_workspace(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Sumc::delete_gsl_sum_levin_u_workspace($self);
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


############# Class : Math::GSL::Sum::gsl_sum_levin_utrunc_workspace ##############

package Math::GSL::Sum::gsl_sum_levin_utrunc_workspace;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Sum );
%OWNER = ();
%ITERATORS = ();
*swig_size_get = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_size_get;
*swig_size_set = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_size_set;
*swig_i_get = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_i_get;
*swig_i_set = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_i_set;
*swig_terms_used_get = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_terms_used_get;
*swig_terms_used_set = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_terms_used_set;
*swig_sum_plain_get = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_sum_plain_get;
*swig_sum_plain_set = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_sum_plain_set;
*swig_q_num_get = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_q_num_get;
*swig_q_num_set = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_q_num_set;
*swig_q_den_get = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_q_den_get;
*swig_q_den_set = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_q_den_set;
*swig_dsum_get = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_dsum_get;
*swig_dsum_set = *Math::GSL::Sumc::gsl_sum_levin_utrunc_workspace_dsum_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Sumc::new_gsl_sum_levin_utrunc_workspace(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Sumc::delete_gsl_sum_levin_utrunc_workspace($self);
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

package Math::GSL::Sum;


@EXPORT_OK = qw/
               gsl_sum_levin_u_alloc 
               gsl_sum_levin_u_free 
               gsl_sum_levin_u_accel 
               gsl_sum_levin_u_minmax 
               gsl_sum_levin_u_step 
               gsl_sum_levin_utrunc_alloc 
               gsl_sum_levin_utrunc_free 
               gsl_sum_levin_utrunc_accel 
               gsl_sum_levin_utrunc_minmax 
               gsl_sum_levin_utrunc_step 
             /;
%EXPORT_TAGS = ( all => \@EXPORT_OK );
1;

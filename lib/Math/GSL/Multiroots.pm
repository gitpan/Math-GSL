# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::Multiroots;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Multirootsc;
bootstrap Math::GSL::Multiroots;
package Math::GSL::Multiroots;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Multiroots;

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

package Math::GSL::Multiroots;

*gsl_multiroot_fdjacobian = *Math::GSL::Multirootsc::gsl_multiroot_fdjacobian;
*gsl_multiroot_fsolver_alloc = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_alloc;
*gsl_multiroot_fsolver_free = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_free;
*gsl_multiroot_fsolver_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_set;
*gsl_multiroot_fsolver_iterate = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_iterate;
*gsl_multiroot_fsolver_name = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_name;
*gsl_multiroot_fsolver_root = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_root;
*gsl_multiroot_fsolver_dx = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_dx;
*gsl_multiroot_fsolver_f = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_f;
*gsl_multiroot_fdfsolver_alloc = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_alloc;
*gsl_multiroot_fdfsolver_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_set;
*gsl_multiroot_fdfsolver_iterate = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_iterate;
*gsl_multiroot_fdfsolver_free = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_free;
*gsl_multiroot_fdfsolver_name = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_name;
*gsl_multiroot_fdfsolver_root = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_root;
*gsl_multiroot_fdfsolver_dx = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_dx;
*gsl_multiroot_fdfsolver_f = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_f;
*gsl_multiroot_test_delta = *Math::GSL::Multirootsc::gsl_multiroot_test_delta;
*gsl_multiroot_test_residual = *Math::GSL::Multirootsc::gsl_multiroot_test_residual;

############# Class : Math::GSL::Multiroots::gsl_multiroot_function_struct ##############

package Math::GSL::Multiroots::gsl_multiroot_function_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multiroots );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Multirootsc::gsl_multiroot_function_struct_f_get;
*swig_f_set = *Math::GSL::Multirootsc::gsl_multiroot_function_struct_f_set;
*swig_n_get = *Math::GSL::Multirootsc::gsl_multiroot_function_struct_n_get;
*swig_n_set = *Math::GSL::Multirootsc::gsl_multiroot_function_struct_n_set;
*swig_params_get = *Math::GSL::Multirootsc::gsl_multiroot_function_struct_params_get;
*swig_params_set = *Math::GSL::Multirootsc::gsl_multiroot_function_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multirootsc::new_gsl_multiroot_function_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multirootsc::delete_gsl_multiroot_function_struct($self);
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


############# Class : Math::GSL::Multiroots::gsl_multiroot_fsolver_type ##############

package Math::GSL::Multiroots::gsl_multiroot_fsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multiroots );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_name_get;
*swig_name_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_name_set;
*swig_size_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_size_get;
*swig_size_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_size_set;
*swig_alloc_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_alloc_get;
*swig_alloc_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_alloc_set;
*swig_set_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_set_get;
*swig_set_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_iterate_set;
*swig_free_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_free_get;
*swig_free_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_free_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multirootsc::new_gsl_multiroot_fsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multirootsc::delete_gsl_multiroot_fsolver_type($self);
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


############# Class : Math::GSL::Multiroots::gsl_multiroot_fsolver ##############

package Math::GSL::Multiroots::gsl_multiroot_fsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multiroots );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_get;
*swig_type_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_type_set;
*swig_function_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_function_get;
*swig_function_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_function_set;
*swig_x_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_x_get;
*swig_x_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_x_set;
*swig_f_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_f_get;
*swig_f_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_f_set;
*swig_dx_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_dx_get;
*swig_dx_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_dx_set;
*swig_state_get = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_state_get;
*swig_state_set = *Math::GSL::Multirootsc::gsl_multiroot_fsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multirootsc::new_gsl_multiroot_fsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multirootsc::delete_gsl_multiroot_fsolver($self);
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


############# Class : Math::GSL::Multiroots::gsl_multiroot_function_fdf_struct ##############

package Math::GSL::Multiroots::gsl_multiroot_function_fdf_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multiroots );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_f_get;
*swig_f_set = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_f_set;
*swig_df_get = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_df_get;
*swig_df_set = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_df_set;
*swig_fdf_get = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_fdf_get;
*swig_fdf_set = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_fdf_set;
*swig_n_get = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_n_get;
*swig_n_set = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_n_set;
*swig_params_get = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_params_get;
*swig_params_set = *Math::GSL::Multirootsc::gsl_multiroot_function_fdf_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multirootsc::new_gsl_multiroot_function_fdf_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multirootsc::delete_gsl_multiroot_function_fdf_struct($self);
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


############# Class : Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type ##############

package Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multiroots );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_name_get;
*swig_name_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_name_set;
*swig_size_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_size_get;
*swig_size_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_size_set;
*swig_alloc_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_alloc_get;
*swig_alloc_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_alloc_set;
*swig_set_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_set_get;
*swig_set_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_iterate_set;
*swig_free_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_free_get;
*swig_free_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_free_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multirootsc::new_gsl_multiroot_fdfsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multirootsc::delete_gsl_multiroot_fdfsolver_type($self);
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


############# Class : Math::GSL::Multiroots::gsl_multiroot_fdfsolver ##############

package Math::GSL::Multiroots::gsl_multiroot_fdfsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Multiroots );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_get;
*swig_type_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_type_set;
*swig_fdf_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_fdf_get;
*swig_fdf_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_fdf_set;
*swig_x_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_x_get;
*swig_x_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_x_set;
*swig_f_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_f_get;
*swig_f_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_f_set;
*swig_J_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_J_get;
*swig_J_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_J_set;
*swig_dx_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_dx_get;
*swig_dx_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_dx_set;
*swig_state_get = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_state_get;
*swig_state_set = *Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Multirootsc::new_gsl_multiroot_fdfsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Multirootsc::delete_gsl_multiroot_fdfsolver($self);
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

package Math::GSL::Multiroots;

*GSL_MAJOR_VERSION = *Math::GSL::Multirootsc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Multirootsc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Multirootsc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Multirootsc::GSL_NEGZERO;

my %__gsl_multiroot_fsolver_dnewton_hash;
tie %__gsl_multiroot_fsolver_dnewton_hash,"Math::GSL::Multiroots::gsl_multiroot_fsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fsolver_dnewton;
$gsl_multiroot_fsolver_dnewton= \%__gsl_multiroot_fsolver_dnewton_hash;
bless $gsl_multiroot_fsolver_dnewton, Math::GSL::Multiroots::gsl_multiroot_fsolver_type;

my %__gsl_multiroot_fsolver_broyden_hash;
tie %__gsl_multiroot_fsolver_broyden_hash,"Math::GSL::Multiroots::gsl_multiroot_fsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fsolver_broyden;
$gsl_multiroot_fsolver_broyden= \%__gsl_multiroot_fsolver_broyden_hash;
bless $gsl_multiroot_fsolver_broyden, Math::GSL::Multiroots::gsl_multiroot_fsolver_type;

my %__gsl_multiroot_fsolver_hybrid_hash;
tie %__gsl_multiroot_fsolver_hybrid_hash,"Math::GSL::Multiroots::gsl_multiroot_fsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fsolver_hybrid;
$gsl_multiroot_fsolver_hybrid= \%__gsl_multiroot_fsolver_hybrid_hash;
bless $gsl_multiroot_fsolver_hybrid, Math::GSL::Multiroots::gsl_multiroot_fsolver_type;

my %__gsl_multiroot_fsolver_hybrids_hash;
tie %__gsl_multiroot_fsolver_hybrids_hash,"Math::GSL::Multiroots::gsl_multiroot_fsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fsolver_hybrids;
$gsl_multiroot_fsolver_hybrids= \%__gsl_multiroot_fsolver_hybrids_hash;
bless $gsl_multiroot_fsolver_hybrids, Math::GSL::Multiroots::gsl_multiroot_fsolver_type;

my %__gsl_multiroot_fdfsolver_newton_hash;
tie %__gsl_multiroot_fdfsolver_newton_hash,"Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_newton;
$gsl_multiroot_fdfsolver_newton= \%__gsl_multiroot_fdfsolver_newton_hash;
bless $gsl_multiroot_fdfsolver_newton, Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type;

my %__gsl_multiroot_fdfsolver_gnewton_hash;
tie %__gsl_multiroot_fdfsolver_gnewton_hash,"Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_gnewton;
$gsl_multiroot_fdfsolver_gnewton= \%__gsl_multiroot_fdfsolver_gnewton_hash;
bless $gsl_multiroot_fdfsolver_gnewton, Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type;

my %__gsl_multiroot_fdfsolver_hybridj_hash;
tie %__gsl_multiroot_fdfsolver_hybridj_hash,"Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_hybridj;
$gsl_multiroot_fdfsolver_hybridj= \%__gsl_multiroot_fdfsolver_hybridj_hash;
bless $gsl_multiroot_fdfsolver_hybridj, Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type;

my %__gsl_multiroot_fdfsolver_hybridsj_hash;
tie %__gsl_multiroot_fdfsolver_hybridsj_hash,"Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type", $Math::GSL::Multirootsc::gsl_multiroot_fdfsolver_hybridsj;
$gsl_multiroot_fdfsolver_hybridsj= \%__gsl_multiroot_fdfsolver_hybridsj_hash;
bless $gsl_multiroot_fdfsolver_hybridsj, Math::GSL::Multiroots::gsl_multiroot_fdfsolver_type;

@EXPORT_OK = qw/
               gsl_multiroot_fdjacobian 
               gsl_multiroot_fsolver_alloc 
               gsl_multiroot_fsolver_free 
               gsl_multiroot_fsolver_set 
               gsl_multiroot_fsolver_iterate 
               gsl_multiroot_fsolver_name 
               gsl_multiroot_fsolver_root 
               gsl_multiroot_fsolver_dx 
               gsl_multiroot_fsolver_f 
               gsl_multiroot_fdfsolver_alloc 
               gsl_multiroot_fdfsolver_set 
               gsl_multiroot_fdfsolver_iterate 
               gsl_multiroot_fdfsolver_free 
               gsl_multiroot_fdfsolver_name 
               gsl_multiroot_fdfsolver_root 
               gsl_multiroot_fdfsolver_dx 
               gsl_multiroot_fdfsolver_f 
               gsl_multiroot_test_delta 
               gsl_multiroot_test_residual 
               $gsl_multiroot_fsolver_dnewton
               $gsl_multiroot_fsolver_broyden
               $gsl_multiroot_fsolver_hybrid
               $gsl_multiroot_fsolver_hybrids
               $gsl_multiroot_fdfsolver_newton
               $gsl_multiroot_fdfsolver_gnewton
               $gsl_multiroot_fdfsolver_hybridj
               $gsl_multiroot_fdfsolver_hybridsj
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Multiroots - Multidimensional root-finding 


=head1 SYNOPSIS

This module is not yet implemented. Patches Welcome!

    use Math::GSL::Multiroots qw/:all/;

Solving nonlinear systems with n equations in n unknowns.

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over 

=item * C<gsl_multiroot_fdjacobian >

=item * C<gsl_multiroot_fsolver_alloc >

=item * C<gsl_multiroot_fsolver_free >

=item * C<gsl_multiroot_fsolver_set >

=item * C<gsl_multiroot_fsolver_iterate >

=item * C<gsl_multiroot_fsolver_name >

=item * C<gsl_multiroot_fsolver_root >

=item * C<gsl_multiroot_fsolver_dx >

=item * C<gsl_multiroot_fsolver_f >

=item * C<gsl_multiroot_fdfsolver_alloc >

=item * C<gsl_multiroot_fdfsolver_set >

=item * C<gsl_multiroot_fdfsolver_iterate >

=item * C<gsl_multiroot_fdfsolver_free >

=item * C<gsl_multiroot_fdfsolver_name >

=item * C<gsl_multiroot_fdfsolver_root >

=item * C<gsl_multiroot_fdfsolver_dx >

=item * C<gsl_multiroot_fdfsolver_f >

=item * C<gsl_multiroot_test_delta >

=item * C<gsl_multiroot_test_residual >

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

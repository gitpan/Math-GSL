# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.31
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Roots;
require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
package Math::GSL::Rootsc;
bootstrap Math::GSL::Roots;
package Math::GSL::Roots;
@EXPORT = qw( );

# ---------- BASE METHODS -------------

package Math::GSL::Roots;

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

package Math::GSL::Roots;

*gsl_root_fsolver_alloc = *Math::GSL::Rootsc::gsl_root_fsolver_alloc;
*gsl_root_fsolver_free = *Math::GSL::Rootsc::gsl_root_fsolver_free;
*gsl_root_fsolver_set = *Math::GSL::Rootsc::gsl_root_fsolver_set;
*gsl_root_fsolver_iterate = *Math::GSL::Rootsc::gsl_root_fsolver_iterate;
*gsl_root_fsolver_name = *Math::GSL::Rootsc::gsl_root_fsolver_name;
*gsl_root_fsolver_root = *Math::GSL::Rootsc::gsl_root_fsolver_root;
*gsl_root_fsolver_x_lower = *Math::GSL::Rootsc::gsl_root_fsolver_x_lower;
*gsl_root_fsolver_x_upper = *Math::GSL::Rootsc::gsl_root_fsolver_x_upper;
*gsl_root_fdfsolver_alloc = *Math::GSL::Rootsc::gsl_root_fdfsolver_alloc;
*gsl_root_fdfsolver_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_set;
*gsl_root_fdfsolver_iterate = *Math::GSL::Rootsc::gsl_root_fdfsolver_iterate;
*gsl_root_fdfsolver_free = *Math::GSL::Rootsc::gsl_root_fdfsolver_free;
*gsl_root_fdfsolver_name = *Math::GSL::Rootsc::gsl_root_fdfsolver_name;
*gsl_root_fdfsolver_root = *Math::GSL::Rootsc::gsl_root_fdfsolver_root;
*gsl_root_test_interval = *Math::GSL::Rootsc::gsl_root_test_interval;
*gsl_root_test_residual = *Math::GSL::Rootsc::gsl_root_test_residual;
*gsl_root_test_delta = *Math::GSL::Rootsc::gsl_root_test_delta;

############# Class : Math::GSL::Roots::gsl_root_fsolver_type ##############

package Math::GSL::Roots::gsl_root_fsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_name_get;
*swig_name_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_name_set;
*swig_size_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_size_get;
*swig_size_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_size_set;
*swig_set_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_set_get;
*swig_set_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_iterate_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fsolver_type($self);
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


############# Class : Math::GSL::Roots::gsl_root_fsolver ##############

package Math::GSL::Roots::gsl_root_fsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Rootsc::gsl_root_fsolver_type_get;
*swig_type_set = *Math::GSL::Rootsc::gsl_root_fsolver_type_set;
*swig_function_get = *Math::GSL::Rootsc::gsl_root_fsolver_function_get;
*swig_function_set = *Math::GSL::Rootsc::gsl_root_fsolver_function_set;
*swig_root_get = *Math::GSL::Rootsc::gsl_root_fsolver_root_get;
*swig_root_set = *Math::GSL::Rootsc::gsl_root_fsolver_root_set;
*swig_x_lower_get = *Math::GSL::Rootsc::gsl_root_fsolver_x_lower_get;
*swig_x_lower_set = *Math::GSL::Rootsc::gsl_root_fsolver_x_lower_set;
*swig_x_upper_get = *Math::GSL::Rootsc::gsl_root_fsolver_x_upper_get;
*swig_x_upper_set = *Math::GSL::Rootsc::gsl_root_fsolver_x_upper_set;
*swig_state_get = *Math::GSL::Rootsc::gsl_root_fsolver_state_get;
*swig_state_set = *Math::GSL::Rootsc::gsl_root_fsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fsolver($self);
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


############# Class : Math::GSL::Roots::gsl_root_fdfsolver_type ##############

package Math::GSL::Roots::gsl_root_fdfsolver_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_name_get;
*swig_name_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_name_set;
*swig_size_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_size_get;
*swig_size_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_size_set;
*swig_set_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_set_get;
*swig_set_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_set_set;
*swig_iterate_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_iterate_get;
*swig_iterate_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_iterate_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fdfsolver_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fdfsolver_type($self);
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


############# Class : Math::GSL::Roots::gsl_root_fdfsolver ##############

package Math::GSL::Roots::gsl_root_fdfsolver;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Roots );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_get;
*swig_type_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_type_set;
*swig_fdf_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_fdf_get;
*swig_fdf_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_fdf_set;
*swig_root_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_root_get;
*swig_root_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_root_set;
*swig_state_get = *Math::GSL::Rootsc::gsl_root_fdfsolver_state_get;
*swig_state_set = *Math::GSL::Rootsc::gsl_root_fdfsolver_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Rootsc::new_gsl_root_fdfsolver(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Rootsc::delete_gsl_root_fdfsolver($self);
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

package Math::GSL::Roots;


my %__gsl_root_fsolver_bisection_hash;
tie %__gsl_root_fsolver_bisection_hash,"Math::GSL::Roots::gsl_root_fsolver_type", $Math::GSL::Rootsc::gsl_root_fsolver_bisection;
$gsl_root_fsolver_bisection= \%__gsl_root_fsolver_bisection_hash;
bless $gsl_root_fsolver_bisection, Math::GSL::Roots::gsl_root_fsolver_type;

my %__gsl_root_fsolver_brent_hash;
tie %__gsl_root_fsolver_brent_hash,"Math::GSL::Roots::gsl_root_fsolver_type", $Math::GSL::Rootsc::gsl_root_fsolver_brent;
$gsl_root_fsolver_brent= \%__gsl_root_fsolver_brent_hash;
bless $gsl_root_fsolver_brent, Math::GSL::Roots::gsl_root_fsolver_type;

my %__gsl_root_fsolver_falsepos_hash;
tie %__gsl_root_fsolver_falsepos_hash,"Math::GSL::Roots::gsl_root_fsolver_type", $Math::GSL::Rootsc::gsl_root_fsolver_falsepos;
$gsl_root_fsolver_falsepos= \%__gsl_root_fsolver_falsepos_hash;
bless $gsl_root_fsolver_falsepos, Math::GSL::Roots::gsl_root_fsolver_type;

my %__gsl_root_fdfsolver_newton_hash;
tie %__gsl_root_fdfsolver_newton_hash,"Math::GSL::Roots::gsl_root_fdfsolver_type", $Math::GSL::Rootsc::gsl_root_fdfsolver_newton;
$gsl_root_fdfsolver_newton= \%__gsl_root_fdfsolver_newton_hash;
bless $gsl_root_fdfsolver_newton, Math::GSL::Roots::gsl_root_fdfsolver_type;

my %__gsl_root_fdfsolver_secant_hash;
tie %__gsl_root_fdfsolver_secant_hash,"Math::GSL::Roots::gsl_root_fdfsolver_type", $Math::GSL::Rootsc::gsl_root_fdfsolver_secant;
$gsl_root_fdfsolver_secant= \%__gsl_root_fdfsolver_secant_hash;
bless $gsl_root_fdfsolver_secant, Math::GSL::Roots::gsl_root_fdfsolver_type;

my %__gsl_root_fdfsolver_steffenson_hash;
tie %__gsl_root_fdfsolver_steffenson_hash,"Math::GSL::Roots::gsl_root_fdfsolver_type", $Math::GSL::Rootsc::gsl_root_fdfsolver_steffenson;
$gsl_root_fdfsolver_steffenson= \%__gsl_root_fdfsolver_steffenson_hash;
bless $gsl_root_fdfsolver_steffenson, Math::GSL::Roots::gsl_root_fdfsolver_type;

@EXPORT_OK = qw/
               gsl_root_fsolver_alloc 
               gsl_root_fsolver_free 
               gsl_root_fsolver_set 
               gsl_root_fsolver_iterate 
               gsl_root_fsolver_name 
               gsl_root_fsolver_root 
               gsl_root_fsolver_x_lower 
               gsl_root_fsolver_x_upper 
               gsl_root_fdfsolver_alloc 
               gsl_root_fdfsolver_set 
               gsl_root_fdfsolver_iterate 
               gsl_root_fdfsolver_free 
               gsl_root_fdfsolver_name 
               gsl_root_fdfsolver_root 
               gsl_root_test_interval 
               gsl_root_test_residual 
               gsl_root_test_delta 
               $gsl_root_fsolver_bisection    
               $gsl_root_fsolver_brent   
               $gsl_root_fsolver_falsepos     
               $gsl_root_fdfsolver_newton     
               $gsl_root_fdfsolver_secant     
               $gsl_root_fdfsolver_steffenson 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Roots - Routines for finding roots of arbitrary one-dimensional functions.

=head1 SYNOPSIS

This module is not yet implemented. Patches Welcome!

use Math::GSL::Roots qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over 

=item * C<gsl_root_fsolver_alloc >

=item * C<gsl_root_fsolver_free >

=item * C<gsl_root_fsolver_set >

=item * C<gsl_root_fsolver_iterate >

=item * C<gsl_root_fsolver_name >

=item * C<gsl_root_fsolver_root >

=item * C<gsl_root_fsolver_x_lower >

=item * C<gsl_root_fsolver_x_upper >

=item * C<gsl_root_fdfsolver_alloc >

=item * C<gsl_root_fdfsolver_set >

=item * C<gsl_root_fdfsolver_iterate >

=item * C<gsl_root_fdfsolver_free >

=item * C<gsl_root_fdfsolver_name >

=item * C<gsl_root_fdfsolver_root >

=item * C<gsl_root_test_interval >

=item * C<gsl_root_test_residual >

=item * C<gsl_root_test_delta >

=back

This module also includes the following constants :

=over

=item * C<$gsl_root_fsolver_bisection>

=item * C<$gsl_root_fsolver_brent>

=item * C<$gsl_root_fsolver_falsepos>

=item * C<$gsl_root_fdfsolver_newton>

=item * C<$gsl_root_fdfsolver_secant>

=item * C<$gsl_root_fdfsolver_steffenson>

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

 Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want


=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

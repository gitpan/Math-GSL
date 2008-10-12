# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.31
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Min;
require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
package Math::GSL::Minc;
bootstrap Math::GSL::Min;
package Math::GSL::Min;
@EXPORT = qw( );

# ---------- BASE METHODS -------------

package Math::GSL::Min;

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

package Math::GSL::Min;

*gsl_min_fminimizer_alloc = *Math::GSL::Minc::gsl_min_fminimizer_alloc;
*gsl_min_fminimizer_free = *Math::GSL::Minc::gsl_min_fminimizer_free;
*gsl_min_fminimizer_set = *Math::GSL::Minc::gsl_min_fminimizer_set;
*gsl_min_fminimizer_set_with_values = *Math::GSL::Minc::gsl_min_fminimizer_set_with_values;
*gsl_min_fminimizer_iterate = *Math::GSL::Minc::gsl_min_fminimizer_iterate;
*gsl_min_fminimizer_name = *Math::GSL::Minc::gsl_min_fminimizer_name;
*gsl_min_fminimizer_x_minimum = *Math::GSL::Minc::gsl_min_fminimizer_x_minimum;
*gsl_min_fminimizer_x_lower = *Math::GSL::Minc::gsl_min_fminimizer_x_lower;
*gsl_min_fminimizer_x_upper = *Math::GSL::Minc::gsl_min_fminimizer_x_upper;
*gsl_min_fminimizer_f_minimum = *Math::GSL::Minc::gsl_min_fminimizer_f_minimum;
*gsl_min_fminimizer_f_lower = *Math::GSL::Minc::gsl_min_fminimizer_f_lower;
*gsl_min_fminimizer_f_upper = *Math::GSL::Minc::gsl_min_fminimizer_f_upper;
*gsl_min_fminimizer_minimum = *Math::GSL::Minc::gsl_min_fminimizer_minimum;
*gsl_min_test_interval = *Math::GSL::Minc::gsl_min_test_interval;
*gsl_min_find_bracket = *Math::GSL::Minc::gsl_min_find_bracket;
*gsl_max = *Math::GSL::Minc::gsl_max;
*gsl_min = *Math::GSL::Minc::gsl_min;

############# Class : Math::GSL::Min::gsl_min_fminimizer_type ##############

package Math::GSL::Min::gsl_min_fminimizer_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Min );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::Minc::gsl_min_fminimizer_type_name_get;
*swig_name_set = *Math::GSL::Minc::gsl_min_fminimizer_type_name_set;
*swig_size_get = *Math::GSL::Minc::gsl_min_fminimizer_type_size_get;
*swig_size_set = *Math::GSL::Minc::gsl_min_fminimizer_type_size_set;
*swig_set_get = *Math::GSL::Minc::gsl_min_fminimizer_type_set_get;
*swig_set_set = *Math::GSL::Minc::gsl_min_fminimizer_type_set_set;
*swig_iterate_get = *Math::GSL::Minc::gsl_min_fminimizer_type_iterate_get;
*swig_iterate_set = *Math::GSL::Minc::gsl_min_fminimizer_type_iterate_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Minc::new_gsl_min_fminimizer_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Minc::delete_gsl_min_fminimizer_type($self);
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


############# Class : Math::GSL::Min::gsl_min_fminimizer ##############

package Math::GSL::Min::gsl_min_fminimizer;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Min );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::Minc::gsl_min_fminimizer_type_get;
*swig_type_set = *Math::GSL::Minc::gsl_min_fminimizer_type_set;
*swig_function_get = *Math::GSL::Minc::gsl_min_fminimizer_function_get;
*swig_function_set = *Math::GSL::Minc::gsl_min_fminimizer_function_set;
*swig_x_minimum_get = *Math::GSL::Minc::gsl_min_fminimizer_x_minimum_get;
*swig_x_minimum_set = *Math::GSL::Minc::gsl_min_fminimizer_x_minimum_set;
*swig_x_lower_get = *Math::GSL::Minc::gsl_min_fminimizer_x_lower_get;
*swig_x_lower_set = *Math::GSL::Minc::gsl_min_fminimizer_x_lower_set;
*swig_x_upper_get = *Math::GSL::Minc::gsl_min_fminimizer_x_upper_get;
*swig_x_upper_set = *Math::GSL::Minc::gsl_min_fminimizer_x_upper_set;
*swig_f_minimum_get = *Math::GSL::Minc::gsl_min_fminimizer_f_minimum_get;
*swig_f_minimum_set = *Math::GSL::Minc::gsl_min_fminimizer_f_minimum_set;
*swig_f_lower_get = *Math::GSL::Minc::gsl_min_fminimizer_f_lower_get;
*swig_f_lower_set = *Math::GSL::Minc::gsl_min_fminimizer_f_lower_set;
*swig_f_upper_get = *Math::GSL::Minc::gsl_min_fminimizer_f_upper_get;
*swig_f_upper_set = *Math::GSL::Minc::gsl_min_fminimizer_f_upper_set;
*swig_state_get = *Math::GSL::Minc::gsl_min_fminimizer_state_get;
*swig_state_set = *Math::GSL::Minc::gsl_min_fminimizer_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Minc::new_gsl_min_fminimizer(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Minc::delete_gsl_min_fminimizer($self);
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


############# Class : Math::GSL::Min::gsl_function_struct ##############

package Math::GSL::Min::gsl_function_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Min );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::Minc::gsl_function_struct_function_get;
*swig_function_set = *Math::GSL::Minc::gsl_function_struct_function_set;
*swig_params_get = *Math::GSL::Minc::gsl_function_struct_params_get;
*swig_params_set = *Math::GSL::Minc::gsl_function_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Minc::new_gsl_function_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Minc::delete_gsl_function_struct($self);
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


############# Class : Math::GSL::Min::gsl_function_fdf_struct ##############

package Math::GSL::Min::gsl_function_fdf_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Min );
%OWNER = ();
%ITERATORS = ();
*swig_f_get = *Math::GSL::Minc::gsl_function_fdf_struct_f_get;
*swig_f_set = *Math::GSL::Minc::gsl_function_fdf_struct_f_set;
*swig_df_get = *Math::GSL::Minc::gsl_function_fdf_struct_df_get;
*swig_df_set = *Math::GSL::Minc::gsl_function_fdf_struct_df_set;
*swig_fdf_get = *Math::GSL::Minc::gsl_function_fdf_struct_fdf_get;
*swig_fdf_set = *Math::GSL::Minc::gsl_function_fdf_struct_fdf_set;
*swig_params_get = *Math::GSL::Minc::gsl_function_fdf_struct_params_get;
*swig_params_set = *Math::GSL::Minc::gsl_function_fdf_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Minc::new_gsl_function_fdf_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Minc::delete_gsl_function_fdf_struct($self);
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


############# Class : Math::GSL::Min::gsl_function_vec_struct ##############

package Math::GSL::Min::gsl_function_vec_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::Min );
%OWNER = ();
%ITERATORS = ();
*swig_function_get = *Math::GSL::Minc::gsl_function_vec_struct_function_get;
*swig_function_set = *Math::GSL::Minc::gsl_function_vec_struct_function_set;
*swig_params_get = *Math::GSL::Minc::gsl_function_vec_struct_params_get;
*swig_params_set = *Math::GSL::Minc::gsl_function_vec_struct_params_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::Minc::new_gsl_function_vec_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::Minc::delete_gsl_function_vec_struct($self);
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

package Math::GSL::Min;


my %__gsl_min_fminimizer_goldensection_hash;
tie %__gsl_min_fminimizer_goldensection_hash,"Math::GSL::Min::gsl_min_fminimizer_type", $Math::GSL::Minc::gsl_min_fminimizer_goldensection;
$gsl_min_fminimizer_goldensection= \%__gsl_min_fminimizer_goldensection_hash;
bless $gsl_min_fminimizer_goldensection, Math::GSL::Min::gsl_min_fminimizer_type;

my %__gsl_min_fminimizer_brent_hash;
tie %__gsl_min_fminimizer_brent_hash,"Math::GSL::Min::gsl_min_fminimizer_type", $Math::GSL::Minc::gsl_min_fminimizer_brent;
$gsl_min_fminimizer_brent= \%__gsl_min_fminimizer_brent_hash;
bless $gsl_min_fminimizer_brent, Math::GSL::Min::gsl_min_fminimizer_type;
*M_E = *Math::GSL::Minc::M_E;
*M_LOG2E = *Math::GSL::Minc::M_LOG2E;
*M_LOG10E = *Math::GSL::Minc::M_LOG10E;
*M_SQRT2 = *Math::GSL::Minc::M_SQRT2;
*M_SQRT1_2 = *Math::GSL::Minc::M_SQRT1_2;
*M_SQRT3 = *Math::GSL::Minc::M_SQRT3;
*M_PI = *Math::GSL::Minc::M_PI;
*M_PI_2 = *Math::GSL::Minc::M_PI_2;
*M_PI_4 = *Math::GSL::Minc::M_PI_4;
*M_SQRTPI = *Math::GSL::Minc::M_SQRTPI;
*M_2_SQRTPI = *Math::GSL::Minc::M_2_SQRTPI;
*M_1_PI = *Math::GSL::Minc::M_1_PI;
*M_2_PI = *Math::GSL::Minc::M_2_PI;
*M_LN10 = *Math::GSL::Minc::M_LN10;
*M_LN2 = *Math::GSL::Minc::M_LN2;
*M_LNPI = *Math::GSL::Minc::M_LNPI;
*M_EULER = *Math::GSL::Minc::M_EULER;


@EXPORT_OK = qw/
   gsl_min_fminimizer_alloc 
   gsl_min_fminimizer_free 
   gsl_min_fminimizer_set 
   gsl_min_fminimizer_set_with_values
   gsl_min_fminimizer_iterate 
   gsl_min_fminimizer_name 
   gsl_min_fminimizer_x_minimum
   gsl_min_fminimizer_x_lower 
   gsl_min_fminimizer_x_upper 
   gsl_min_fminimizer_f_minimum
   gsl_min_fminimizer_f_lower 
   gsl_min_fminimizer_f_upper 
   gsl_min_fminimizer_minimum 
   gsl_min_test_interval 
   gsl_min_find_bracket 
   $gsl_min_fminimizer_brent
   $gsl_min_fminimizer_goldensection
/;

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Min - One-dimensional Minimization

=head1 SYNOPSIS

This module is not yet implemented. Patches Welcome!

use Math::GSL::Min qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over

=item * C<gsl_min_fminimizer_alloc >

=item * C<gsl_min_fminimizer_free >

=item * C<gsl_min_fminimizer_set >

=item * C<gsl_min_fminimizer_set_with_values>

=item * C<gsl_min_fminimizer_iterate >

=item * C<gsl_min_fminimizer_name >

=item * C<gsl_min_fminimizer_x_minimum>

=item * C<gsl_min_fminimizer_x_lower >

=item * C<gsl_min_fminimizer_x_upper >

=item * C<gsl_min_fminimizer_f_minimum>

=item * C<gsl_min_fminimizer_f_lower >

=item * C<gsl_min_fminimizer_f_upper >

=item * C<gsl_min_fminimizer_minimum >

=item * C<gsl_min_test_interval >

=item * C<gsl_min_find_bracket >

=back

This module also includes the following constants :

=over

=item * C<$gsl_min_fminimizer_brent>

=item * C<$gsl_min_fminimizer_goldensection> 

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

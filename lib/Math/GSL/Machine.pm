# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Machine;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Machinec;
bootstrap Math::GSL::Machine;
package Math::GSL::Machine;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Machine;

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

package Math::GSL::Machine;


# ------- VARIABLE STUBS --------

package Math::GSL::Machine;

*GSL_DBL_EPSILON = *Math::GSL::Machinec::GSL_DBL_EPSILON;
*GSL_SQRT_DBL_EPSILON = *Math::GSL::Machinec::GSL_SQRT_DBL_EPSILON;
*GSL_ROOT3_DBL_EPSILON = *Math::GSL::Machinec::GSL_ROOT3_DBL_EPSILON;
*GSL_ROOT4_DBL_EPSILON = *Math::GSL::Machinec::GSL_ROOT4_DBL_EPSILON;
*GSL_ROOT5_DBL_EPSILON = *Math::GSL::Machinec::GSL_ROOT5_DBL_EPSILON;
*GSL_ROOT6_DBL_EPSILON = *Math::GSL::Machinec::GSL_ROOT6_DBL_EPSILON;
*GSL_LOG_DBL_EPSILON = *Math::GSL::Machinec::GSL_LOG_DBL_EPSILON;
*GSL_DBL_MIN = *Math::GSL::Machinec::GSL_DBL_MIN;
*GSL_SQRT_DBL_MIN = *Math::GSL::Machinec::GSL_SQRT_DBL_MIN;
*GSL_ROOT3_DBL_MIN = *Math::GSL::Machinec::GSL_ROOT3_DBL_MIN;
*GSL_ROOT4_DBL_MIN = *Math::GSL::Machinec::GSL_ROOT4_DBL_MIN;
*GSL_ROOT5_DBL_MIN = *Math::GSL::Machinec::GSL_ROOT5_DBL_MIN;
*GSL_ROOT6_DBL_MIN = *Math::GSL::Machinec::GSL_ROOT6_DBL_MIN;
*GSL_LOG_DBL_MIN = *Math::GSL::Machinec::GSL_LOG_DBL_MIN;
*GSL_DBL_MAX = *Math::GSL::Machinec::GSL_DBL_MAX;
*GSL_SQRT_DBL_MAX = *Math::GSL::Machinec::GSL_SQRT_DBL_MAX;
*GSL_ROOT3_DBL_MAX = *Math::GSL::Machinec::GSL_ROOT3_DBL_MAX;
*GSL_ROOT4_DBL_MAX = *Math::GSL::Machinec::GSL_ROOT4_DBL_MAX;
*GSL_ROOT5_DBL_MAX = *Math::GSL::Machinec::GSL_ROOT5_DBL_MAX;
*GSL_ROOT6_DBL_MAX = *Math::GSL::Machinec::GSL_ROOT6_DBL_MAX;
*GSL_LOG_DBL_MAX = *Math::GSL::Machinec::GSL_LOG_DBL_MAX;
*GSL_FLT_EPSILON = *Math::GSL::Machinec::GSL_FLT_EPSILON;
*GSL_SQRT_FLT_EPSILON = *Math::GSL::Machinec::GSL_SQRT_FLT_EPSILON;
*GSL_ROOT3_FLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT3_FLT_EPSILON;
*GSL_ROOT4_FLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT4_FLT_EPSILON;
*GSL_ROOT5_FLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT5_FLT_EPSILON;
*GSL_ROOT6_FLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT6_FLT_EPSILON;
*GSL_LOG_FLT_EPSILON = *Math::GSL::Machinec::GSL_LOG_FLT_EPSILON;
*GSL_FLT_MIN = *Math::GSL::Machinec::GSL_FLT_MIN;
*GSL_SQRT_FLT_MIN = *Math::GSL::Machinec::GSL_SQRT_FLT_MIN;
*GSL_ROOT3_FLT_MIN = *Math::GSL::Machinec::GSL_ROOT3_FLT_MIN;
*GSL_ROOT4_FLT_MIN = *Math::GSL::Machinec::GSL_ROOT4_FLT_MIN;
*GSL_ROOT5_FLT_MIN = *Math::GSL::Machinec::GSL_ROOT5_FLT_MIN;
*GSL_ROOT6_FLT_MIN = *Math::GSL::Machinec::GSL_ROOT6_FLT_MIN;
*GSL_LOG_FLT_MIN = *Math::GSL::Machinec::GSL_LOG_FLT_MIN;
*GSL_FLT_MAX = *Math::GSL::Machinec::GSL_FLT_MAX;
*GSL_SQRT_FLT_MAX = *Math::GSL::Machinec::GSL_SQRT_FLT_MAX;
*GSL_ROOT3_FLT_MAX = *Math::GSL::Machinec::GSL_ROOT3_FLT_MAX;
*GSL_ROOT4_FLT_MAX = *Math::GSL::Machinec::GSL_ROOT4_FLT_MAX;
*GSL_ROOT5_FLT_MAX = *Math::GSL::Machinec::GSL_ROOT5_FLT_MAX;
*GSL_ROOT6_FLT_MAX = *Math::GSL::Machinec::GSL_ROOT6_FLT_MAX;
*GSL_LOG_FLT_MAX = *Math::GSL::Machinec::GSL_LOG_FLT_MAX;
*GSL_SFLT_EPSILON = *Math::GSL::Machinec::GSL_SFLT_EPSILON;
*GSL_SQRT_SFLT_EPSILON = *Math::GSL::Machinec::GSL_SQRT_SFLT_EPSILON;
*GSL_ROOT3_SFLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT3_SFLT_EPSILON;
*GSL_ROOT4_SFLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT4_SFLT_EPSILON;
*GSL_ROOT5_SFLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT5_SFLT_EPSILON;
*GSL_ROOT6_SFLT_EPSILON = *Math::GSL::Machinec::GSL_ROOT6_SFLT_EPSILON;
*GSL_LOG_SFLT_EPSILON = *Math::GSL::Machinec::GSL_LOG_SFLT_EPSILON;
*GSL_MACH_EPS = *Math::GSL::Machinec::GSL_MACH_EPS;
*GSL_SQRT_MACH_EPS = *Math::GSL::Machinec::GSL_SQRT_MACH_EPS;
*GSL_ROOT3_MACH_EPS = *Math::GSL::Machinec::GSL_ROOT3_MACH_EPS;
*GSL_ROOT4_MACH_EPS = *Math::GSL::Machinec::GSL_ROOT4_MACH_EPS;
*GSL_ROOT5_MACH_EPS = *Math::GSL::Machinec::GSL_ROOT5_MACH_EPS;
*GSL_ROOT6_MACH_EPS = *Math::GSL::Machinec::GSL_ROOT6_MACH_EPS;
*GSL_LOG_MACH_EPS = *Math::GSL::Machinec::GSL_LOG_MACH_EPS;

@EXPORT_OK = qw/
               $GSL_DBL_EPSILON 
               $GSL_SQRT_DBL_EPSILON 
               $GSL_ROOT3_DBL_EPSILON 
               $GSL_ROOT4_DBL_EPSILON 
               $GSL_ROOT5_DBL_EPSILON 
               $GSL_ROOT6_DBL_EPSILON 
               $GSL_LOG_DBL_EPSILON 
               $GSL_DBL_MIN 
               $GSL_SQRT_DBL_MIN 
               $GSL_ROOT3_DBL_MIN 
               $GSL_ROOT4_DBL_MIN 
               $GSL_ROOT5_DBL_MIN 
               $GSL_ROOT6_DBL_MIN 
               $GSL_LOG_DBL_MIN 
               $GSL_DBL_MAX 
               $GSL_SQRT_DBL_MAX 
               $GSL_ROOT3_DBL_MAX 
               $GSL_ROOT4_DBL_MAX 
               $GSL_ROOT5_DBL_MAX 
               $GSL_ROOT6_DBL_MAX 
               $GSL_LOG_DBL_MAX 
               $GSL_FLT_EPSILON 
               $GSL_SQRT_FLT_EPSILON 
               $GSL_ROOT3_FLT_EPSILON 
               $GSL_ROOT4_FLT_EPSILON 
               $GSL_ROOT5_FLT_EPSILON 
               $GSL_ROOT6_FLT_EPSILON 
               $GSL_LOG_FLT_EPSILON 
               $GSL_FLT_MIN 
               $GSL_SQRT_FLT_MIN 
               $GSL_ROOT3_FLT_MIN 
               $GSL_ROOT4_FLT_MIN 
               $GSL_ROOT5_FLT_MIN 
               $GSL_ROOT6_FLT_MIN 
               $GSL_LOG_FLT_MIN 
               $GSL_FLT_MAX 
               $GSL_SQRT_FLT_MAX 
               $GSL_ROOT3_FLT_MAX 
               $GSL_ROOT4_FLT_MAX 
               $GSL_ROOT5_FLT_MAX 
               $GSL_ROOT6_FLT_MAX 
               $GSL_LOG_FLT_MAX 
               $GSL_SFLT_EPSILON 
               $GSL_SQRT_SFLT_EPSILON 
               $GSL_ROOT3_SFLT_EPSILON 
               $GSL_ROOT4_SFLT_EPSILON 
               $GSL_ROOT5_SFLT_EPSILON 
               $GSL_ROOT6_SFLT_EPSILON 
               $GSL_LOG_SFLT_EPSILON 
               $GSL_MACH_EPS 
               $GSL_SQRT_MACH_EPS 
               $GSL_ROOT3_MACH_EPS 
               $GSL_ROOT4_MACH_EPS 
               $GSL_ROOT5_MACH_EPS 
               $GSL_ROOT6_MACH_EPS 
               $GSL_LOG_MACH_EPS 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );
__END__

=head1 NAME

Math::GSL::Machine - Machine constants 

=head1 SYNOPSIS

    use Math::GSL::Machine qw/:all/;
    print "The GSL_DBL_EPSILON Constant is $GSL_DBL_EPSILON\n";

Here is a list of all the constants included in this module :

=over

=item *   $GSL_DBL_EPSILON 

=item *   $GSL_SQRT_DBL_EPSILON 

=item *   $GSL_ROOT3_DBL_EPSILON 

=item *   $GSL_ROOT4_DBL_EPSILON 

=item *   $GSL_ROOT5_DBL_EPSILON 

=item *   $GSL_ROOT6_DBL_EPSILON 

=item *   $GSL_LOG_DBL_EPSILON 

=item *   $GSL_DBL_MIN 

=item *   $GSL_SQRT_DBL_MIN 

=item *   $GSL_ROOT3_DBL_MIN 

=item *   $GSL_ROOT4_DBL_MIN 

=item *   $GSL_ROOT5_DBL_MIN 

=item *   $GSL_ROOT6_DBL_MIN 

=item *   $GSL_LOG_DBL_MIN 

=item *   $GSL_DBL_MAX 

=item *   $GSL_SQRT_DBL_MAX 

=item *   $GSL_ROOT3_DBL_MAX 

=item *   $GSL_ROOT4_DBL_MAX 

=item *   $GSL_ROOT5_DBL_MAX 

=item *   $GSL_ROOT6_DBL_MAX 

=item *   $GSL_LOG_DBL_MAX 

=item *   $GSL_FLT_EPSILON 

=item *   $GSL_SQRT_FLT_EPSILON 

=item *   $GSL_ROOT3_FLT_EPSILON 

=item *   $GSL_ROOT4_FLT_EPSILON 

=item *   $GSL_ROOT5_FLT_EPSILON 

=item *   $GSL_ROOT6_FLT_EPSILON 

=item *   $GSL_LOG_FLT_EPSILON 

=item *   $GSL_FLT_MIN 

=item *   $GSL_SQRT_FLT_MIN 

=item *   $GSL_ROOT3_FLT_MIN 

=item *   $GSL_ROOT4_FLT_MIN 

=item *   $GSL_ROOT5_FLT_MIN 

=item *   $GSL_ROOT6_FLT_MIN 

=item *   $GSL_LOG_FLT_MIN 

=item *   $GSL_FLT_MAX 

=item *   $GSL_SQRT_FLT_MAX 

=item *   $GSL_ROOT3_FLT_MAX 

=item *   $GSL_ROOT4_FLT_MAX 

=item *   $GSL_ROOT5_FLT_MAX 

=item *   $GSL_ROOT6_FLT_MAX 

=item *   $GSL_LOG_FLT_MAX 

=item *   $GSL_SFLT_EPSILON 

=item *   $GSL_SQRT_SFLT_EPSILON 

=item *   $GSL_ROOT3_SFLT_EPSILON 

=item *   $GSL_ROOT4_SFLT_EPSILON 

=item *   $GSL_ROOT5_SFLT_EPSILON 

=item *   $GSL_ROOT6_SFLT_EPSILON 

=item *   $GSL_LOG_SFLT_EPSILON 

=item *   $GSL_MACH_EPS 

=item *   $GSL_SQRT_MACH_EPS 

=item *   $GSL_ROOT3_MACH_EPS 

=item *   $GSL_ROOT4_MACH_EPS 

=item *   $GSL_ROOT5_MACH_EPS 

=item *   $GSL_ROOT6_MACH_EPS 

=item *   $GSL_LOG_MACH_EPS 

=back

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

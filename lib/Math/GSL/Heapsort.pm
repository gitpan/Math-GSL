# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.36
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Heapsort;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Heapsortc;
bootstrap Math::GSL::Heapsort;
package Math::GSL::Heapsort;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Heapsort;

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

package Math::GSL::Heapsort;

*gsl_heapsort = *Math::GSL::Heapsortc::gsl_heapsort;
*gsl_heapsort_index = *Math::GSL::Heapsortc::gsl_heapsort_index;

# ------- VARIABLE STUBS --------

package Math::GSL::Heapsort;


@EXPORT_OK = qw/
               gsl_heapsort 
               gsl_heapsort_index 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Heapsort - Functions for sorting data, both directly and indirectly (using an index)

=head1 SYNOPSIS

This module is not yet implemented. Patches Welcome!

use Math::GSL::Heapsort qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over

=item * gsl_heapsort 

=item * gsl_heapsort_index 

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

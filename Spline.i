%module "Math::GSL::Spline"

%include "typemaps.i"
%apply double *OUTPUT { double * y, double * d, double * d2, double * result };

%{
    #include "gsl/gsl_spline.h"
%}

%include "gsl/gsl_spline.h"


%perlcode %{
@EXPORT_OK = qw/
               gsl_spline_alloc 
               gsl_spline_init 
               gsl_spline_name 
               gsl_spline_min_size 
               gsl_spline_eval_e 
               gsl_spline_eval 
               gsl_spline_eval_deriv_e 
               gsl_spline_eval_deriv 
               gsl_spline_eval_deriv2_e 
               gsl_spline_eval_deriv2 
               gsl_spline_eval_integ_e 
               gsl_spline_eval_integ 
               gsl_spline_free 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Spline - Higher-level Interface to Interp

=head1 SYNOPSIS

use Math::GSL::Spline qw /:all/;

=head1 DESCRIPTION

The functions described in the Interp module required the user to supply pointers to the x and y arrays on each call. The following functions are equivalent to the corresponding gsl_interp functions but maintain a copy of this data in the gsl_spline object. This removes the need to pass both xa and ya as arguments on each evaluation.

Here is a list of all the functions in this module :

=over

=item * C<gsl_spline_alloc($T, $size)>

=item * C<gsl_spline_init($spline, $xa, $ya, $size)>

=item * C<gsl_spline_free($spline)>

=item * C<gsl_spline_name($spline)>

=item * C<gsl_spline_min_size($spline)>

=item * C<gsl_spline_eval_e($spline, $x, $acc)>

=item * C<gsl_spline_eval($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv_e($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv2_e($spline, $x, $acc)>

=item * C<gsl_spline_eval_deriv2($spline, $x, $acc)>

=item * C<gsl_spline_eval_integ_e($spline, $a, $b, $acc)>

=item * C<gsl_spline_eval_integ($spline, $a, $b, $acc)>

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


%}

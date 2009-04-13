# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Fit;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Fitc;
bootstrap Math::GSL::Fit;
package Math::GSL::Fit;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Fit;

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

package Math::GSL::Fit;

*gsl_fit_linear = *Math::GSL::Fitc::gsl_fit_linear;
*gsl_fit_wlinear = *Math::GSL::Fitc::gsl_fit_wlinear;
*gsl_fit_linear_est = *Math::GSL::Fitc::gsl_fit_linear_est;
*gsl_fit_mul = *Math::GSL::Fitc::gsl_fit_mul;
*gsl_fit_wmul = *Math::GSL::Fitc::gsl_fit_wmul;
*gsl_fit_mul_est = *Math::GSL::Fitc::gsl_fit_mul_est;

# ------- VARIABLE STUBS --------

package Math::GSL::Fit;

*GSL_MAJOR_VERSION = *Math::GSL::Fitc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::Fitc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::Fitc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::Fitc::GSL_NEGZERO;

@EXPORT_OK = qw/
               gsl_fit_linear 
               gsl_fit_wlinear 
               gsl_fit_linear_est 
               gsl_fit_mul 
               gsl_fit_wmul 
               gsl_fit_mul_est 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::Fit - Least-squares functions for a general linear model with one- or two-parameter regression

=head1 SYNOPSIS

use Math::GSL::Fit qw /:all/;

=head1 DESCRIPTION

The functions in this module perform least-squares fits to a general linear model, y = X c where y is a vector of n observations, X is an n by p matrix of predictor variables, and the elements of the vector c are the p unknown best-fit parameters which are to be estimated.

Here is a list of all the functions in this module :

=over 

=item C<gsl_fit_linear($x, $xstride, $y, $ystride, $n)> - This function computes the best-fit linear regression coefficients (c0,c1) of the model Y = c_0 + c_1 X for the dataset ($x, $y), two vectors (in form of arrays) of length $n with strides $xstride and $ystride. The errors on y are assumed unknown so the variance-covariance matrix for the parameters (c0, c1) is estimated from the scatter of the points around the best-fit line and returned via the parameters (cov00, cov01, cov11). The sum of squares of the residuals from the best-fit line is returned in sumsq. Note: the correlation coefficient of the data can be computed using gsl_stats_correlation (see Correlation), it does not depend on the fit. The function returns the following values in this order : 0 if the operation succeeded, 1 otherwise, c0, c1, cov00, cov01, cov11 and sumsq.

=item C<gsl_fit_wlinear($x, $xstride, $w, $wstride, $y, $ystride, $n)> - This function computes the best-fit linear regression coefficients (c0,c1) of the model Y = c_0 + c_1 X for the weighted dataset ($x, $y), two vectors (in form of arrays) of length $n with strides $xstride and $ystride. The vector (also in the form of an array) $w, of length $n and stride $wstride, specifies the weight of each datapoint. The weight is the reciprocal of the variance for each datapoint in y. The covariance matrix for the parameters (c0, c1) is computed using the weights and returned via the parameters (cov00, cov01, cov11). The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq. The function returns the following values in this order : 0 if the operation succeeded, 1 otherwise, c0, c1, cov00, cov01, cov11 and sumsq.

=item C<gsl_fit_linear_est($x, $c0, $c1, $cov00, $cov01, $cov11)> - This function uses the best-fit linear regression coefficients $c0, $c1 and their covariance $cov00, $cov01, $cov11 to compute the fitted function y and its standard deviation y_err for the model Y = c_0 + c_1 X at the point $x. The function returns the following values in this order : 0 if the operation succeeded, 1 otherwise, y and y_err.

=item C<gsl_fit_mul($x, $xstride, $y, $ystride, $n)> - This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for the datasets ($x, $y), two vectors (in form of arrays) of length $n with strides $xstride and $ystride. The errors on y are assumed unknown so the variance of the parameter c1 is estimated from the scatter of the points around the best-fit line and returned via the parameter cov11. The sum of squares of the residuals from the best-fit line is returned in sumsq. The function returns the following values in this order : 0 if the operation succeeded, 1 otherwise, c1, cov11 and sumsq.

=item C<gsl_fit_wmul($x, $xstride, $w, $wstride, $y, $ystride, $n)> - This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for the weighted datasets ($x, $y), two vectors (in form of arrays) of length $n with strides $xstride and $ystride. The vector (also in the form of an array) $w, of length $n and stride $wstride, specifies the weight of each datapoint. The weight is the reciprocal of the variance for each datapoint in y. The variance of the parameter c1 is computed using the weights and returned via the parameter cov11. The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq. The function returns the following values in this order : 0 if the operation succeeded, 1 otherwise, c1, cov11 and sumsq.

=item C<gsl_fit_mul_est($x, $c1, $cov11)> - This function uses the best-fit linear regression coefficient $c1 and its covariance $cov11 to compute the fitted function y and its standard deviation y_err for the model Y = c_1 X at the point $x. The function returns the following values in this order : 0 if the operation succeeded, 1 otherwise, y and y_err.

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want

=head1 EXAMPLES

This example shows how to use the function gsl_fit_linear. It's important to see that the array passed to to function must be an array reference, not a simple array. Also when you use strides, you need to initialize all the value in the range used, unless you will get warnings.

 my @norris_x = (0.2, 337.4, 118.2, 884.6, 10.1, 226.5, 666.3, 996.3,
                        448.6, 777.0, 558.2, 0.4, 0.6, 775.5, 666.9, 338.0, 
                        447.5, 11.6, 556.0, 228.1, 995.8, 887.6, 120.2, 0.3, 
                        0.3, 556.8, 339.1, 887.2, 999.0, 779.0, 11.1, 118.3,
                        229.2, 669.1, 448.9, 0.5 ) ;
 my @norris_y = ( 0.1, 338.8, 118.1, 888.0, 9.2, 228.1, 668.5, 998.5,
                        449.1, 778.9, 559.2, 0.3, 0.1, 778.1, 668.8, 339.3, 
                        448.9, 10.8, 557.7, 228.3, 998.0, 888.8, 119.6, 0.3, 
                        0.6, 557.6, 339.3, 888.0, 998.5, 778.9, 10.2, 117.6,
                        228.9, 668.4, 449.2, 0.2);
    my $xstride = 2; 
    my $wstride = 3; 
    my $ystride = 5;
    my ($x, $w, $y);
    for my $i (0 .. 175)
    {
        $x->[$i] = 0;
        $w->[$i] = 0;
        $y->[$i] = 0;
    }

    for my $i (0 .. 35)
    {
        $x->[$i*$xstride] = $norris_x[$i];
        $w->[$i*$wstride] = 1.0;
        $y->[$i*$ystride] = $norris_y[$i];
    }
    my ($status, @results) = gsl_fit_linear($x, $xstride, $y, $ystride, 36); 

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

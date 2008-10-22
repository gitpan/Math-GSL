%module "Math::GSL::Statistics"

%include "typemaps.i"
%include "gsl_typemaps.i"

%apply double *OUTPUT { double * min, double * max };

%apply int *OUTPUT { size_t * min_index, size_t * max_index };

%{
    #include "gsl/gsl_statistics_double.h"
    #include "gsl/gsl_statistics_int.h"
    #include "gsl/gsl_statistics_char.h"
%}

%include "gsl/gsl_statistics_double.h"
%include "gsl/gsl_statistics_int.h"
%include "gsl/gsl_statistics_char.h"


%perlcode %{
@EXPORT_OK = qw/
               gsl_stats_mean 
               gsl_stats_variance 
               gsl_stats_sd 
               gsl_stats_variance_with_fixed_mean 
               gsl_stats_sd_with_fixed_mean 
               gsl_stats_tss 
               gsl_stats_tss_m 
               gsl_stats_absdev 
               gsl_stats_skew 
               gsl_stats_kurtosis 
               gsl_stats_lag1_autocorrelation 
               gsl_stats_covariance 
               gsl_stats_correlation 
               gsl_stats_variance_m 
               gsl_stats_sd_m 
               gsl_stats_absdev_m 
               gsl_stats_skew_m_sd 
               gsl_stats_kurtosis_m_sd 
               gsl_stats_lag1_autocorrelation_m 
               gsl_stats_covariance_m 
               gsl_stats_wmean 
               gsl_stats_wvariance 
               gsl_stats_wsd 
               gsl_stats_wvariance_with_fixed_mean 
               gsl_stats_wsd_with_fixed_mean 
               gsl_stats_wtss 
               gsl_stats_wtss_m 
               gsl_stats_wabsdev 
               gsl_stats_wskew 
               gsl_stats_wkurtosis 
               gsl_stats_wvariance_m 
               gsl_stats_wsd_m 
               gsl_stats_wabsdev_m 
               gsl_stats_wskew_m_sd 
               gsl_stats_wkurtosis_m_sd 
               gsl_stats_pvariance 
               gsl_stats_ttest 
               gsl_stats_max 
               gsl_stats_min 
               gsl_stats_minmax 
               gsl_stats_max_index 
               gsl_stats_min_index 
               gsl_stats_minmax_index 
               gsl_stats_median_from_sorted_data 
               gsl_stats_quantile_from_sorted_data 
               /;
our @EXPORT_int = qw/
               gsl_stats_int_mean 
               gsl_stats_int_variance 
               gsl_stats_int_sd 
               gsl_stats_int_variance_with_fixed_mean 
               gsl_stats_int_sd_with_fixed_mean 
               gsl_stats_int_tss 
               gsl_stats_int_tss_m 
               gsl_stats_int_absdev 
               gsl_stats_int_skew 
               gsl_stats_int_kurtosis 
               gsl_stats_int_lag1_autocorrelation 
               gsl_stats_int_covariance 
               gsl_stats_int_correlation 
               gsl_stats_int_variance_m 
               gsl_stats_int_sd_m 
               gsl_stats_int_absdev_m 
               gsl_stats_int_skew_m_sd 
               gsl_stats_int_kurtosis_m_sd 
               gsl_stats_int_lag1_autocorrelation_m 
               gsl_stats_int_covariance_m 
               gsl_stats_int_pvariance 
               gsl_stats_int_ttest 
               gsl_stats_int_max 
               gsl_stats_int_min 
               gsl_stats_int_minmax 
               gsl_stats_int_max_index 
               gsl_stats_int_min_index 
               gsl_stats_int_minmax_index 
               gsl_stats_int_median_from_sorted_data 
               gsl_stats_int_quantile_from_sorted_data 
               /;
our @EXPORT_char = qw/
               gsl_stats_char_mean 
               gsl_stats_char_variance 
               gsl_stats_char_sd 
               gsl_stats_char_variance_with_fixed_mean 
               gsl_stats_char_sd_with_fixed_mean 
               gsl_stats_char_tss 
               gsl_stats_char_tss_m 
               gsl_stats_char_absdev 
               gsl_stats_char_skew 
               gsl_stats_char_kurtosis 
               gsl_stats_char_lag1_autocorrelation 
               gsl_stats_char_covariance 
               gsl_stats_char_correlation 
               gsl_stats_char_variance_m 
               gsl_stats_char_sd_m 
               gsl_stats_char_absdev_m 
               gsl_stats_char_skew_m_sd 
               gsl_stats_char_kurtosis_m_sd 
               gsl_stats_char_lag1_autocorrelation_m 
               gsl_stats_char_covariance_m 
               gsl_stats_char_pvariance 
               gsl_stats_char_ttest 
               gsl_stats_char_max 
               gsl_stats_char_min 
               gsl_stats_char_minmax 
               gsl_stats_char_max_index 
               gsl_stats_char_min_index 
               gsl_stats_char_minmax_index 
               gsl_stats_char_median_from_sorted_data 
               gsl_stats_char_quantile_from_sorted_data 
             /;
push @EXPORT_OK, @EXPORT_int, @EXPORT_char;
%EXPORT_TAGS = ( 
                all => \@EXPORT_OK,
                int => \@EXPORT_int,
                char => \@EXPORT_char
               );
__END__

=head1 NAME

Math::GSL::Statistics - Statistical functions 

=head1 SYNOPSIS

    use Math::GSL::Statistics qw /:all/;
        
    my $data     = [17.2, 18.1, 16.5, 18.3, 12.6];
    my $mean     = gsl_stats_mean($data, 1, 5);
    my $variance = gsl_stats_variance($data, 1, 5);
    my $largest  = gsl_stats_max($data, 1, 5);
    my $smallest = gsl_stats_min($data, 1, 5);
    print qq{
    Dataset : @$data
    Sample mean           $mean 
    Estimated variance    $variance
    Largest value         $largest
    Smallest value        $smallest
    };


=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over 2

=item * C<gsl_stats_mean($data, $stride, $n)> - This function returns the arithmetic mean of the array reference $data, a dataset of length $n with stride $stride. The arithmetic mean, or sample mean, is denoted by \Hat\mu and defined as, \Hat\mu = (1/N) \sum x_i where x_i are the elements of the dataset $data. For samples drawn from a gaussian distribution the variance of \Hat\mu is \sigma^2 / N. 

=item * C<gsl_stats_variance($data, $stride, $n)> - This function returns the estimated, or sample, variance of data, an array reference of length $n with stride $stride. The estimated variance is denoted by \Hat\sigma^2 and is defined by, \Hat\sigma^2 = (1/(N-1)) \sum (x_i - \Hat\mu)^2 where x_i are the elements of the dataset data. Note that the normalization factor of 1/(N-1) results from the derivation of \Hat\sigma^2 as an unbiased estimator of the population variance \sigma^2. For samples drawn from a gaussian distribution the variance of \Hat\sigma^2 itself is 2 \sigma^4 / N. This function computes the mean via a call to gsl_stats_mean. If you have already computed the mean then you can pass it directly to gsl_stats_variance_m. 

=item * C<gsl_stats_sd($data, $stride, $n)>

=item * C<gsl_stats_sd_m($data, $stride, $n, $mean)>

The standard deviation is defined as the square root of the variance. These functions return the square root of the corresponding variance functions above.

=item * C<gsl_stats_variance_with_fixed_mean($data, $stride, $n, $mean)> - This function calculates the standard deviation of the array reference $data for a fixed population mean $mean. The result is the square root of the corresponding variance function.

=item * C<gsl_stats_sd_with_fixed_mean($data, $stride, $n, $mean)> - This function computes an unbiased estimate of the variance of data when the population mean $mean of the underlying distribution is known a priori. In this case the estimator for the variance uses the factor 1/N and the sample mean \Hat\mu is replaced by the known population mean \mu, \Hat\sigma^2 = (1/N) \sum (x_i - \mu)^2

=item * C<gsl_stats_tss($data, $stride, $n)>

=item * C<gsl_stats_tss_m($data, $stride, $n, $mean)>

These functions return the total sum of squares (TSS) of data about the mean. For gsl_stats_tss_m the user-supplied value of mean is used, and for gsl_stats_tss it is computed using gsl_stats_mean. TSS =  \sum (x_i - mean)^2

=item * C<gsl_stats_absdev($data, $stride, $n)> - This function computes the absolute deviation from the mean of data, a dataset of length $n with stride $stride. The absolute deviation from the mean is defined as, absdev  = (1/N) \sum |x_i - \Hat\mu| where x_i are the elements of the array reference $data. The absolute deviation from the mean provides a more robust measure of the width of a distribution than the variance. This function computes the mean of data via a call to gsl_stats_mean. 

=item * C<gsl_stats_skew($data, $stride, $n)> - This function computes the skewness of $data, a dataset in the form of an array reference of length $n with stride $stride. The skewness is defined as, skew = (1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^3 where x_i are the elements of the dataset $data. The skewness measures the asymmetry of the tails of a distribution. The function computes the mean and estimated standard deviation of data via calls to gsl_stats_mean and gsl_stats_sd. 

=item * C<gsl_stats_skew_m_sd($data, $stride, $n, $mean, $sd)> - This function computes the skewness of the array reference $data using the given values of the mean $mean and standard deviation $sd, skew = (1/N) \sum ((x_i - mean)/sd)^3. These functions are useful if you have already computed the mean and standard deviation of $data and want to avoid recomputing them. 

=item * C<gsl_stats_kurtosis($data, $stride, $n)> - This function computes the kurtosis of data, an array reference of length $n with stride $stride. The kurtosis is defined as, kurtosis = ((1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^4)  - 3. The kurtosis measures how sharply peaked a distribution is, relative to its width. The kurtosis is normalized to zero for a gaussian distribution. 

=item * C<gsl_stats_kurtosis_m_sd($data, $stride, $n, $mean, $sd)> - This function computes the kurtosis of the array reference $data using the given values of the mean $mean and standard deviation $sd, kurtosis = ((1/N) \sum ((x_i - mean)/sd)^4) - 3. This function is useful if you have already computed the mean and standard deviation of data and want to avoid recomputing them. 

=item * C<gsl_stats_lag1_autocorrelation($data, $stride, $n)> - This function computes the lag-1 autocorrelation of the array reference data.
 a_1 = {\sum_{i = 1}^{n} (x_{i} - \Hat\mu) (x_{i-1} - \Hat\mu)
  \over
 \sum_{i = 1}^{n} (x_{i} - \Hat\mu) (x_{i} - \Hat\mu)}

=item * C<gsl_stats_lag1_autocorrelation_m($data, $stride, $n, $mean)> - This function computes the lag-1 autocorrelation of the array reference $data using the given value of the mean $mean.

=item * C<gsl_stats_covariance($data1, $stride1, $data2, $stride2, $n)> - This function computes the covariance of the array reference $data1 and $data2 which must both be of the same length $n. covar = (1/(n - 1)) \sum_{i = 1}^{n} (x_i - \Hat x) (y_i - \Hat y)

=item * C<gsl_stats_covariance_m($data1, $stride1, $data2, $stride2, $n, $mean1, $mean2)> - This function computes the covariance of the array reference $data1 and $data2 using the given values of the means, $mean1 and $mean2. This is useful if you have already computed the means of $data1 and $data2 and want to avoid recomputing them.

=item * C<gsl_stats_correlation($data1, $stride1, $data2, $stride2, $n)> - This function efficiently computes the Pearson correlation coefficient between the array reference $data1 and $data2 which must both be of the same length $n.
 r = cov(x, y) / (\Hat\sigma_x \Hat\sigma_y)
   = {1/(n-1) \sum (x_i - \Hat x) (y_i - \Hat y)
      \over
      \sqrt{1/(n-1) \sum (x_i - \Hat x)^2} \sqrt{1/(n-1) \sum (y_i - \Hat y)^2}
     }
         
=item * C<gsl_stats_variance_m($data, $stride, $n, $mean)> - This function returns the sample variance of $data, an array reference, relative to the given value of $mean. The function is computed with \Hat\mu replaced by the value of mean that you supply, \Hat\sigma^2 = (1/(N-1)) \sum (x_i - mean)^2

=item * C<gsl_stats_absdev_m($data, $stride, $n, $mean)> - This function computes the absolute deviation of the dataset $data, an array refrence, relative to the given value of $mean, absdev  = (1/N) \sum |x_i - mean|. This function is useful if you have already computed the mean of data (and want to avoid recomputing it), or wish to calculate the absolute deviation relative to another value (such as zero, or the median). 

=item * C<gsl_stats_wmean($w, $wstride, $data, $stride, $n)> - This function returns the weighted mean of the dataset $data array reference with stride $stride and length $n, using the set of weights $w, which is an array reference, with stride $wstride and length $n. The weighted mean is defined as, \Hat\mu = (\sum w_i x_i) / (\sum w_i)

=item * C<gsl_stats_wvariance($w, $wstride, $data, $stride, $n)> - This function returns the estimated variance of the dataset $data, which is the dataset, with stride $stride and length $n, using the set of weights $w (as an array reference) with stride $wstride and length $n. The estimated variance of a weighted dataset is defined as,  \Hat\sigma^2 = ((\sum w_i)/((\sum w_i)^2 - \sum (w_i^2))) \sum w_i (x_i - \Hat\mu)^2. Note that this expression reduces to an unweighted variance with the familiar 1/(N-1) factor when there are N equal non-zero weights. 

=item * C<gsl_stats_wvariance_m($w, $wstride, $data, $stride, $n, $wmean, $wsd)> - This function returns the estimated variance of the weighted dataset $data (which is an array reference) using the given weighted mean $wmean.

=item * C<gsl_stats_wsd($w, $wstride, $data, $stride, $n)> - The standard deviation is defined as the square root of the variance. This function returns the square root of the corresponding variance function gsl_stats_wvariance above.

=item * C<gsl_stats_wsd_m($w, $wstride, $data, $stride, $n, $wmean)> - This function returns the square root of the corresponding variance function gsl_stats_wvariance_m above.

=item * C<gsl_stats_wvariance_with_fixed_mean($w, $wstride, $data, $stride, $n, $mean)> - This function computes an unbiased estimate of the variance of weighted dataset $data (which is an array reference) when the population mean $mean of the underlying distribution is known a priori. In this case the estimator for the variance replaces the sample mean \Hat\mu by the known population mean \mu, \Hat\sigma^2 = (\sum w_i (x_i - \mu)^2) / (\sum w_i)

=item * C<gsl_stats_wsd_with_fixed_mean($w, $wstride, $data, $stride, $n, $mean)> - The standard deviation is defined as the square root of the variance. This function returns the square root of the corresponding variance function above.

=item * C<gsl_stats_wtss($w, $wstride, $data, $stride, $n)>

=item * C<gsl_stats_wtss_m($w, $wstride, $data, $stride, $n, $wmean)> - These functions return the weighted total sum of squares (TSS) of data about the weighted mean. For gsl_stats_wtss_m the user-supplied value of $wmean is used, and for gsl_stats_wtss it is computed using gsl_stats_wmean. TSS =  \sum w_i (x_i - wmean)^2

=item * C<gsl_stats_wabsdev($w, $wstride, $data, $stride, $n)> - This function computes the weighted absolute deviation from the weighted mean of $data, which is an array reference. The absolute deviation from the mean is defined as, absdev = (\sum w_i |x_i - \Hat\mu|) / (\sum w_i)

=item * C<gsl_stats_wabsdev_m($w, $wstride, $data, $stride, $n, $wmean)> - This function computes the absolute deviation of the weighted dataset $data (an array reference) about the given weighted mean $wmean.

=item * C<gsl_stats_wskew($w, $wstride, $data, $stride, $n)> - This function computes the weighted skewness of the dataset $data, an array reference. skew = (\sum w_i ((x_i - xbar)/\sigma)^3) / (\sum w_i)

=item * C<gsl_stats_wskew_m_sd($w, $wstride, $data, $stride, $n, $wmean, $wsd)> - This function computes the weighted skewness of the dataset $data using the given values of the weighted mean and weighted standard deviation, $wmean and $wsd.

=item * C<gsl_stats_wkurtosis($w, $wstride, $data, $stride, $n)> - This function computes the weighted kurtosis of the dataset $data, an array reference. kurtosis = ((\sum w_i ((x_i - xbar)/sigma)^4) / (\sum w_i)) - 3

=item * C<gsl_stats_wkurtosis_m_sd($w, $wstride, $data, $stride, $n, $wmean, $wsd)> - This function computes the weighted kurtosis of the dataset $data, an array reference, using the given values of the weighted mean and weighted standard deviation, $wmean and $wsd. 

=item * C<gsl_stats_pvariance($data, $stride, $n, $data2, $stride2, $n2)>

=item * C<gsl_stats_ttest($data1, $stride1, $n1, $data2, $stride2, $n2)>

=item * C<gsl_stats_max($data, $stride, $n)> - This function returns the maximum value in the $data array reference, a dataset of length $n with stride $stride. The maximum value is defined as the value of the element x_i which satisfies x_i >= x_j for all j. If you want instead to find the element with the largest absolute magnitude you will need to apply fabs or abs to your data before calling this function. 

=item * C<gsl_stats_min($data, $stride, $n)> - This function returns the minimum value in $data (which is an array reference) a dataset of length $n with stride $stride. The minimum value is defined as the value of the element x_i which satisfies x_i <= x_j for all j. If you want instead to find the element with the smallest absolute magnitude you will need to apply fabs or abs to your data before calling this function. 

=item * C<gsl_stats_minmax($data, $stride, $n)> - This function finds both the minimum and maximum values in $data, which is an array reference, in a single pass and returns them in this order.

=item * C<gsl_stats_max_index($data, $stride, $n)> - This function returns the index of the maximum value in $data array reference, a dataset of length $n with stride $stride. The maximum value is defined as the value of the element x_i which satisfies x_i >= x_j for all j. When there are several equal maximum elements then the first one is chosen.

=item * C<gsl_stats_min_index($data, $stride, $n)> - This function returns the index of the minimum value in $data array reference, a dataset of length $n with stride $stride. The minimum value is defined as the value of the element x_i which satisfies x_i <= x_j for all j. When there are several equal minimum elements then the first one is chosen.

=item * C<gsl_stats_minmax_index($data, $stride, $n)> - This function returns the indexes of the minimum and maximum values in $data, an array reference in a single pass. The value are returned in this order. 

=item * C<gsl_stats_median_from_sorted_data($sorted_data, $stride, $n)> - This function returns the median value of $sorted_data (which is an array reference), a dataset of length $n with stride $stride. The elements of the array must be in ascending numerical order. There are no checks to see whether the data are sorted, so the function gsl_sort should always be used first. This function can be found in the Math::GSL::Sort module.  When the dataset has an odd number of elements the median is the value of element (n-1)/2. When the dataset has an even number of elements the median is the mean of the two nearest middle values, elements (n-1)/2 and n/2. Since the algorithm for computing the median involves interpolation this function always returns a floating-point number, even for integer data types.

=item * C<gsl_stats_quantile_from_sorted_data($sorted_data, $stride, $n, $f)> - This function returns a quantile value of $sorted_data, a double-precision array reference of length $n with stride $stride. The elements of the array must be in ascending numerical order. The quantile is determined by the f, a fraction between 0 and 1. For example, to compute the value of the 75th percentile f should have the value 0.75. There are no checks to see whether the data are sorted, so the function gsl_sort should always be used first. This function can be found in the Math::GSL::Sort module. The quantile is found by interpolation, using the formula quantile = (1 - \delta) x_i + \delta x_{i+1} where i is floor((n - 1)f) and \delta is (n-1)f - i. Thus the minimum value of the array (data[0*stride]) is given by f equal to zero, the maximum value (data[(n-1)*stride]) is given by f equal to one and the median value is given by f equal to 0.5. Since the algorithm for computing quantiles involves interpolation this function always returns a floating-point number, even for integer data types. 

=back

The following function are simply variants for int and char of the last functions:

=over 4

=item * C<gsl_stats_int_mean >

=item * C<gsl_stats_int_variance >

=item * C<gsl_stats_int_sd >

=item * C<gsl_stats_int_variance_with_fixed_mean >

=item * C<gsl_stats_int_sd_with_fixed_mean >

=item * C<gsl_stats_int_tss >

=item * C<gsl_stats_int_tss_m >

=item * C<gsl_stats_int_absdev >

=item * C<gsl_stats_int_skew >

=item * C<gsl_stats_int_kurtosis >

=item * C<gsl_stats_int_lag1_autocorrelation >

=item * C<gsl_stats_int_covariance >

=item * C<gsl_stats_int_correlation >

=item * C<gsl_stats_int_variance_m >

=item * C<gsl_stats_int_sd_m >

=item * C<gsl_stats_int_absdev_m >

=item * C<gsl_stats_int_skew_m_sd >

=item * C<gsl_stats_int_kurtosis_m_sd >

=item * C<gsl_stats_int_lag1_autocorrelation_m >

=item * C<gsl_stats_int_covariance_m >

=item * C<gsl_stats_int_pvariance >

=item * C<gsl_stats_int_ttest >

=item * C<gsl_stats_int_max >

=item * C<gsl_stats_int_min >

=item * C<gsl_stats_int_minmax >

=item * C<gsl_stats_int_max_index >

=item * C<gsl_stats_int_min_index >

=item * C<gsl_stats_int_minmax_index >

=item * C<gsl_stats_int_median_from_sorted_data >

=item * C<gsl_stats_int_quantile_from_sorted_data >

=item * C<gsl_stats_char_mean >

=item * C<gsl_stats_char_variance >

=item * C<gsl_stats_char_sd >

=item * C<gsl_stats_char_variance_with_fixed_mean >

=item * C<gsl_stats_char_sd_with_fixed_mean >

=item * C<gsl_stats_char_tss >

=item * C<gsl_stats_char_tss_m >

=item * C<gsl_stats_char_absdev >

=item * C<gsl_stats_char_skew >

=item * C<gsl_stats_char_kurtosis >

=item * C<gsl_stats_char_lag1_autocorrelation >

=item * C<gsl_stats_char_covariance >

=item * C<gsl_stats_char_correlation >

=item * C<gsl_stats_char_variance_m >

=item * C<gsl_stats_char_sd_m >

=item * C<gsl_stats_char_absdev_m >

=item * C<gsl_stats_char_skew_m_sd >

=item * C<gsl_stats_char_kurtosis_m_sd >

=item * C<gsl_stats_char_lag1_autocorrelation_m >

=item * C<gsl_stats_char_covariance_m >

=item * C<gsl_stats_char_pvariance >

=item * C<gsl_stats_char_ttest >

=item * C<gsl_stats_char_max >

=item * C<gsl_stats_char_min >

=item * C<gsl_stats_char_minmax >

=item * C<gsl_stats_char_max_index >

=item * C<gsl_stats_char_min_index >

=item * C<gsl_stats_char_minmax_index >

=item * C<gsl_stats_char_median_from_sorted_data >

=item * C<gsl_stats_char_quantile_from_sorted_data >

=back

You have to add the functions you want to use inside the qw /put_funtion_here /. 
You can also write use Math::GSL::Statistics qw/:all/; to use all avaible functions of the module. 
Other tags are also avaible, here is a complete list of all tags for this module :

=over

=item all

=item int

=item char

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

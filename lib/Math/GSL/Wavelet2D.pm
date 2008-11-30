# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::Wavelet2D;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::Wavelet2Dc;
bootstrap Math::GSL::Wavelet2D;
package Math::GSL::Wavelet2D;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::Wavelet2D;

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

package Math::GSL::Wavelet2D;

*gsl_wavelet2d_transform = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_transform;
*gsl_wavelet2d_transform_forward = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_transform_forward;
*gsl_wavelet2d_transform_inverse = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_transform_inverse;
*gsl_wavelet2d_nstransform = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_nstransform;
*gsl_wavelet2d_nstransform_forward = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_nstransform_forward;
*gsl_wavelet2d_nstransform_inverse = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_nstransform_inverse;
*gsl_wavelet2d_transform_matrix = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_transform_matrix;
*gsl_wavelet2d_transform_matrix_forward = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_transform_matrix_forward;
*gsl_wavelet2d_transform_matrix_inverse = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_transform_matrix_inverse;
*gsl_wavelet2d_nstransform_matrix = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_nstransform_matrix;
*gsl_wavelet2d_nstransform_matrix_forward = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_nstransform_matrix_forward;
*gsl_wavelet2d_nstransform_matrix_inverse = *Math::GSL::Wavelet2Dc::gsl_wavelet2d_nstransform_matrix_inverse;

# ------- VARIABLE STUBS --------

package Math::GSL::Wavelet2D;



@EXPORT_OK = qw/
    gsl_wavelet2d_transform 
    gsl_wavelet2d_transform_forward 
    gsl_wavelet2d_transform_inverse 
    gsl_wavelet2d_nstransform 
    gsl_wavelet2d_nstransform_forward 
    gsl_wavelet2d_nstransform_inverse 
    gsl_wavelet2d_transform_matrix 
    gsl_wavelet2d_transform_matrix_forward 
    gsl_wavelet2d_transform_matrix_inverse 
    gsl_wavelet2d_nstransform_matrix 
    gsl_wavelet2d_nstransform_matrix_forward 
    gsl_wavelet2d_nstransform_matrix_inverse 

/;

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

1;

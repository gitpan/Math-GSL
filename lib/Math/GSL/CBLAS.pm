# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::CBLAS;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::CBLASc;
bootstrap Math::GSL::CBLAS;
package Math::GSL::CBLAS;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::CBLAS;

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

package Math::GSL::CBLAS;

*cblas_sdsdot = *Math::GSL::CBLASc::cblas_sdsdot;
*cblas_dsdot = *Math::GSL::CBLASc::cblas_dsdot;
*cblas_sdot = *Math::GSL::CBLASc::cblas_sdot;
*cblas_ddot = *Math::GSL::CBLASc::cblas_ddot;
*cblas_cdotu_sub = *Math::GSL::CBLASc::cblas_cdotu_sub;
*cblas_cdotc_sub = *Math::GSL::CBLASc::cblas_cdotc_sub;
*cblas_zdotu_sub = *Math::GSL::CBLASc::cblas_zdotu_sub;
*cblas_zdotc_sub = *Math::GSL::CBLASc::cblas_zdotc_sub;
*cblas_snrm2 = *Math::GSL::CBLASc::cblas_snrm2;
*cblas_sasum = *Math::GSL::CBLASc::cblas_sasum;
*cblas_dnrm2 = *Math::GSL::CBLASc::cblas_dnrm2;
*cblas_dasum = *Math::GSL::CBLASc::cblas_dasum;
*cblas_scnrm2 = *Math::GSL::CBLASc::cblas_scnrm2;
*cblas_scasum = *Math::GSL::CBLASc::cblas_scasum;
*cblas_dznrm2 = *Math::GSL::CBLASc::cblas_dznrm2;
*cblas_dzasum = *Math::GSL::CBLASc::cblas_dzasum;
*cblas_isamax = *Math::GSL::CBLASc::cblas_isamax;
*cblas_idamax = *Math::GSL::CBLASc::cblas_idamax;
*cblas_icamax = *Math::GSL::CBLASc::cblas_icamax;
*cblas_izamax = *Math::GSL::CBLASc::cblas_izamax;
*cblas_sswap = *Math::GSL::CBLASc::cblas_sswap;
*cblas_scopy = *Math::GSL::CBLASc::cblas_scopy;
*cblas_saxpy = *Math::GSL::CBLASc::cblas_saxpy;
*cblas_dswap = *Math::GSL::CBLASc::cblas_dswap;
*cblas_dcopy = *Math::GSL::CBLASc::cblas_dcopy;
*cblas_daxpy = *Math::GSL::CBLASc::cblas_daxpy;
*cblas_cswap = *Math::GSL::CBLASc::cblas_cswap;
*cblas_ccopy = *Math::GSL::CBLASc::cblas_ccopy;
*cblas_caxpy = *Math::GSL::CBLASc::cblas_caxpy;
*cblas_zswap = *Math::GSL::CBLASc::cblas_zswap;
*cblas_zcopy = *Math::GSL::CBLASc::cblas_zcopy;
*cblas_zaxpy = *Math::GSL::CBLASc::cblas_zaxpy;
*cblas_srotg = *Math::GSL::CBLASc::cblas_srotg;
*cblas_srotmg = *Math::GSL::CBLASc::cblas_srotmg;
*cblas_srot = *Math::GSL::CBLASc::cblas_srot;
*cblas_srotm = *Math::GSL::CBLASc::cblas_srotm;
*cblas_drotg = *Math::GSL::CBLASc::cblas_drotg;
*cblas_drotmg = *Math::GSL::CBLASc::cblas_drotmg;
*cblas_drot = *Math::GSL::CBLASc::cblas_drot;
*cblas_drotm = *Math::GSL::CBLASc::cblas_drotm;
*cblas_sscal = *Math::GSL::CBLASc::cblas_sscal;
*cblas_dscal = *Math::GSL::CBLASc::cblas_dscal;
*cblas_cscal = *Math::GSL::CBLASc::cblas_cscal;
*cblas_zscal = *Math::GSL::CBLASc::cblas_zscal;
*cblas_csscal = *Math::GSL::CBLASc::cblas_csscal;
*cblas_zdscal = *Math::GSL::CBLASc::cblas_zdscal;
*cblas_sgemv = *Math::GSL::CBLASc::cblas_sgemv;
*cblas_sgbmv = *Math::GSL::CBLASc::cblas_sgbmv;
*cblas_strmv = *Math::GSL::CBLASc::cblas_strmv;
*cblas_stbmv = *Math::GSL::CBLASc::cblas_stbmv;
*cblas_stpmv = *Math::GSL::CBLASc::cblas_stpmv;
*cblas_strsv = *Math::GSL::CBLASc::cblas_strsv;
*cblas_stbsv = *Math::GSL::CBLASc::cblas_stbsv;
*cblas_stpsv = *Math::GSL::CBLASc::cblas_stpsv;
*cblas_dgemv = *Math::GSL::CBLASc::cblas_dgemv;
*cblas_dgbmv = *Math::GSL::CBLASc::cblas_dgbmv;
*cblas_dtrmv = *Math::GSL::CBLASc::cblas_dtrmv;
*cblas_dtbmv = *Math::GSL::CBLASc::cblas_dtbmv;
*cblas_dtpmv = *Math::GSL::CBLASc::cblas_dtpmv;
*cblas_dtrsv = *Math::GSL::CBLASc::cblas_dtrsv;
*cblas_dtbsv = *Math::GSL::CBLASc::cblas_dtbsv;
*cblas_dtpsv = *Math::GSL::CBLASc::cblas_dtpsv;
*cblas_cgemv = *Math::GSL::CBLASc::cblas_cgemv;
*cblas_cgbmv = *Math::GSL::CBLASc::cblas_cgbmv;
*cblas_ctrmv = *Math::GSL::CBLASc::cblas_ctrmv;
*cblas_ctbmv = *Math::GSL::CBLASc::cblas_ctbmv;
*cblas_ctpmv = *Math::GSL::CBLASc::cblas_ctpmv;
*cblas_ctrsv = *Math::GSL::CBLASc::cblas_ctrsv;
*cblas_ctbsv = *Math::GSL::CBLASc::cblas_ctbsv;
*cblas_ctpsv = *Math::GSL::CBLASc::cblas_ctpsv;
*cblas_zgemv = *Math::GSL::CBLASc::cblas_zgemv;
*cblas_zgbmv = *Math::GSL::CBLASc::cblas_zgbmv;
*cblas_ztrmv = *Math::GSL::CBLASc::cblas_ztrmv;
*cblas_ztbmv = *Math::GSL::CBLASc::cblas_ztbmv;
*cblas_ztpmv = *Math::GSL::CBLASc::cblas_ztpmv;
*cblas_ztrsv = *Math::GSL::CBLASc::cblas_ztrsv;
*cblas_ztbsv = *Math::GSL::CBLASc::cblas_ztbsv;
*cblas_ztpsv = *Math::GSL::CBLASc::cblas_ztpsv;
*cblas_ssymv = *Math::GSL::CBLASc::cblas_ssymv;
*cblas_ssbmv = *Math::GSL::CBLASc::cblas_ssbmv;
*cblas_sspmv = *Math::GSL::CBLASc::cblas_sspmv;
*cblas_sger = *Math::GSL::CBLASc::cblas_sger;
*cblas_ssyr = *Math::GSL::CBLASc::cblas_ssyr;
*cblas_sspr = *Math::GSL::CBLASc::cblas_sspr;
*cblas_ssyr2 = *Math::GSL::CBLASc::cblas_ssyr2;
*cblas_sspr2 = *Math::GSL::CBLASc::cblas_sspr2;
*cblas_dsymv = *Math::GSL::CBLASc::cblas_dsymv;
*cblas_dsbmv = *Math::GSL::CBLASc::cblas_dsbmv;
*cblas_dspmv = *Math::GSL::CBLASc::cblas_dspmv;
*cblas_dger = *Math::GSL::CBLASc::cblas_dger;
*cblas_dsyr = *Math::GSL::CBLASc::cblas_dsyr;
*cblas_dspr = *Math::GSL::CBLASc::cblas_dspr;
*cblas_dsyr2 = *Math::GSL::CBLASc::cblas_dsyr2;
*cblas_dspr2 = *Math::GSL::CBLASc::cblas_dspr2;
*cblas_chemv = *Math::GSL::CBLASc::cblas_chemv;
*cblas_chbmv = *Math::GSL::CBLASc::cblas_chbmv;
*cblas_chpmv = *Math::GSL::CBLASc::cblas_chpmv;
*cblas_cgeru = *Math::GSL::CBLASc::cblas_cgeru;
*cblas_cgerc = *Math::GSL::CBLASc::cblas_cgerc;
*cblas_cher = *Math::GSL::CBLASc::cblas_cher;
*cblas_chpr = *Math::GSL::CBLASc::cblas_chpr;
*cblas_cher2 = *Math::GSL::CBLASc::cblas_cher2;
*cblas_chpr2 = *Math::GSL::CBLASc::cblas_chpr2;
*cblas_zhemv = *Math::GSL::CBLASc::cblas_zhemv;
*cblas_zhbmv = *Math::GSL::CBLASc::cblas_zhbmv;
*cblas_zhpmv = *Math::GSL::CBLASc::cblas_zhpmv;
*cblas_zgeru = *Math::GSL::CBLASc::cblas_zgeru;
*cblas_zgerc = *Math::GSL::CBLASc::cblas_zgerc;
*cblas_zher = *Math::GSL::CBLASc::cblas_zher;
*cblas_zhpr = *Math::GSL::CBLASc::cblas_zhpr;
*cblas_zher2 = *Math::GSL::CBLASc::cblas_zher2;
*cblas_zhpr2 = *Math::GSL::CBLASc::cblas_zhpr2;
*cblas_sgemm = *Math::GSL::CBLASc::cblas_sgemm;
*cblas_ssymm = *Math::GSL::CBLASc::cblas_ssymm;
*cblas_ssyrk = *Math::GSL::CBLASc::cblas_ssyrk;
*cblas_ssyr2k = *Math::GSL::CBLASc::cblas_ssyr2k;
*cblas_strmm = *Math::GSL::CBLASc::cblas_strmm;
*cblas_strsm = *Math::GSL::CBLASc::cblas_strsm;
*cblas_dgemm = *Math::GSL::CBLASc::cblas_dgemm;
*cblas_dsymm = *Math::GSL::CBLASc::cblas_dsymm;
*cblas_dsyrk = *Math::GSL::CBLASc::cblas_dsyrk;
*cblas_dsyr2k = *Math::GSL::CBLASc::cblas_dsyr2k;
*cblas_dtrmm = *Math::GSL::CBLASc::cblas_dtrmm;
*cblas_dtrsm = *Math::GSL::CBLASc::cblas_dtrsm;
*cblas_cgemm = *Math::GSL::CBLASc::cblas_cgemm;
*cblas_csymm = *Math::GSL::CBLASc::cblas_csymm;
*cblas_csyrk = *Math::GSL::CBLASc::cblas_csyrk;
*cblas_csyr2k = *Math::GSL::CBLASc::cblas_csyr2k;
*cblas_ctrmm = *Math::GSL::CBLASc::cblas_ctrmm;
*cblas_ctrsm = *Math::GSL::CBLASc::cblas_ctrsm;
*cblas_zgemm = *Math::GSL::CBLASc::cblas_zgemm;
*cblas_zsymm = *Math::GSL::CBLASc::cblas_zsymm;
*cblas_zsyrk = *Math::GSL::CBLASc::cblas_zsyrk;
*cblas_zsyr2k = *Math::GSL::CBLASc::cblas_zsyr2k;
*cblas_ztrmm = *Math::GSL::CBLASc::cblas_ztrmm;
*cblas_ztrsm = *Math::GSL::CBLASc::cblas_ztrsm;
*cblas_chemm = *Math::GSL::CBLASc::cblas_chemm;
*cblas_cherk = *Math::GSL::CBLASc::cblas_cherk;
*cblas_cher2k = *Math::GSL::CBLASc::cblas_cher2k;
*cblas_zhemm = *Math::GSL::CBLASc::cblas_zhemm;
*cblas_zherk = *Math::GSL::CBLASc::cblas_zherk;
*cblas_zher2k = *Math::GSL::CBLASc::cblas_zher2k;

# ------- VARIABLE STUBS --------

package Math::GSL::CBLAS;

*GSL_MAJOR_VERSION = *Math::GSL::CBLASc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::CBLASc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::CBLASc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::CBLASc::GSL_NEGZERO;
*CblasRowMajor = *Math::GSL::CBLASc::CblasRowMajor;
*CblasColMajor = *Math::GSL::CBLASc::CblasColMajor;
*CblasNoTrans = *Math::GSL::CBLASc::CblasNoTrans;
*CblasTrans = *Math::GSL::CBLASc::CblasTrans;
*CblasConjTrans = *Math::GSL::CBLASc::CblasConjTrans;
*CblasUpper = *Math::GSL::CBLASc::CblasUpper;
*CblasLower = *Math::GSL::CBLASc::CblasLower;
*CblasNonUnit = *Math::GSL::CBLASc::CblasNonUnit;
*CblasUnit = *Math::GSL::CBLASc::CblasUnit;
*CblasLeft = *Math::GSL::CBLASc::CblasLeft;
*CblasRight = *Math::GSL::CBLASc::CblasRight;

@EXPORT_OK = qw/
               cblas_sdsdot 
               cblas_dsdot 
               cblas_sdot 
               cblas_ddot 
               cblas_cdotu_sub 
               cblas_cdotc_sub 
               cblas_zdotu_sub 
               cblas_zdotc_sub 
               cblas_snrm2 
               cblas_sasum 
               cblas_dnrm2 
               cblas_dasum 
               cblas_scnrm2 
               cblas_scasum 
               cblas_dznrm2 
               cblas_dzasum 
               cblas_isamax 
               cblas_idamax 
               cblas_icamax 
               cblas_izamax 
               cblas_sswap 
               cblas_scopy 
               cblas_saxpy 
               cblas_dswap 
               cblas_dcopy 
               cblas_daxpy 
               cblas_cswap 
               cblas_ccopy 
               cblas_caxpy 
               cblas_zswap 
               cblas_zcopy 
               cblas_zaxpy 
               cblas_srotg 
               cblas_srotmg 
               cblas_srot 
               cblas_srotm 
               cblas_drotg 
               cblas_drotmg 
               cblas_drot 
               cblas_drotm 
               cblas_sscal 
               cblas_dscal 
               cblas_cscal 
               cblas_zscal 
               cblas_csscal 
               cblas_zdscal 
               cblas_sgemv 
               cblas_sgbmv 
               cblas_strmv 
               cblas_stbmv 
               cblas_stpmv 
               cblas_strsv 
               cblas_stbsv 
               cblas_stpsv 
               cblas_dgemv 
               cblas_dgbmv 
               cblas_dtrmv 
               cblas_dtbmv 
               cblas_dtpmv 
               cblas_dtrsv 
               cblas_dtbsv 
               cblas_dtpsv 
               cblas_cgemv 
               cblas_cgbmv 
               cblas_ctrmv 
               cblas_ctbmv 
               cblas_ctpmv 
               cblas_ctrsv 
               cblas_ctbsv 
               cblas_ctpsv 
               cblas_zgemv 
               cblas_zgbmv 
               cblas_ztrmv 
               cblas_ztbmv 
               cblas_ztpmv 
               cblas_ztrsv 
               cblas_ztbsv 
               cblas_ztpsv 
               cblas_ssymv 
               cblas_ssbmv 
               cblas_sspmv 
               cblas_sger 
               cblas_ssyr 
               cblas_sspr 
               cblas_ssyr2 
               cblas_sspr2 
               cblas_dsymv 
               cblas_dsbmv 
               cblas_dspmv 
               cblas_dger 
               cblas_dsyr 
               cblas_dspr 
               cblas_dsyr2 
               cblas_dspr2 
               cblas_chemv 
               cblas_chbmv 
               cblas_chpmv 
               cblas_cgeru 
               cblas_cgerc 
               cblas_cher 
               cblas_chpr 
               cblas_cher2 
               cblas_chpr2 
               cblas_zhemv 
               cblas_zhbmv 
               cblas_zhpmv 
               cblas_zgeru 
               cblas_zgerc 
               cblas_zher 
               cblas_zhpr 
               cblas_zher2 
               cblas_zhpr2 
               cblas_sgemm 
               cblas_ssymm 
               cblas_ssyrk 
               cblas_ssyr2k 
               cblas_strmm 
               cblas_strsm 
               cblas_dgemm 
               cblas_dsymm 
               cblas_dsyrk 
               cblas_dsyr2k 
               cblas_dtrmm 
               cblas_dtrsm 
               cblas_cgemm 
               cblas_csymm 
               cblas_csyrk 
               cblas_csyr2k 
               cblas_ctrmm 
               cblas_ctrsm 
               cblas_zgemm 
               cblas_zsymm 
               cblas_zsyrk 
               cblas_zsyr2k 
               cblas_ztrmm 
               cblas_ztrsm 
               cblas_chemm 
               cblas_cherk 
               cblas_cher2k 
               cblas_zhemm 
               cblas_zherk 
               cblas_zher2k 
               cblas_xerbla 
               $CblasRowMajor 
               $CblasColMajor 
               $CblasNoTrans 
               $CblasTrans 
               $CblasConjTrans 
               $CblasUpper 
               $CblasLower 
               $CblasNonUnit 
               $CblasUnit 
               $CblasLeft 
               $CblasRight 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::CBLAS - Basic Linear Algebra Subprograms based on C functions

=head1 SYNOPSIS

use Math::GSL::CBLAS qw/:all/;

=head1 DESCRIPTION

Here is a list of all the functions included in this module :

=over 1

=item C<cblas_sdsdot>

=item C<cblas_dsdot>

=item C<cblas_sdot>

=item C<cblas_ddot>

=item C<cblas_cdotu_sub>

=item C<cblas_cdotc_sub>

=item C<cblas_zdotu_sub>

=item C<cblas_zdotc_sub>

=item C<cblas_snrm2>

=item C<cblas_sasum>

=item C<cblas_dnrm2>

=item C<cblas_dasum>

=item C<cblas_scnrm2>

=item C<cblas_scasum>

=item C<cblas_dznrm2>

=item C<cblas_dzasum>

=item C<cblas_isamax>

=item C<cblas_idamax>

=item C<cblas_icamax>

=item C<cblas_izamax>

=item C<cblas_sswap>

=item C<cblas_scopy>

=item C<cblas_saxpy>

=item C<cblas_dswap>

=item C<cblas_dcopy>

=item C<cblas_daxpy>

=item C<cblas_cswap>

=item C<cblas_ccopy>

=item C<cblas_caxpy>

=item C<cblas_zswap>

=item C<cblas_zcopy>

=item C<cblas_zaxpy>

=item C<cblas_srotg>

=item C<cblas_srotmg>

=item C<cblas_srot>

=item C<cblas_srotm>

=item C<cblas_drotg>

=item C<cblas_drotmg>

=item C<cblas_drot>

=item C<cblas_drotm>

=item C<cblas_sscal>

=item C<cblas_dscal>

=item C<cblas_cscal>

=item C<cblas_zscal>

=item C<cblas_csscal>

=item C<cblas_zdscal>

=item C<cblas_sgemv>

=item C<cblas_sgbmv>

=item C<cblas_strmv>

=item C<cblas_stbmv>

=item C<cblas_stpmv>

=item C<cblas_strsv>

=item C<cblas_stbsv>

=item C<cblas_stpsv>

=item C<cblas_dgemv>

=item C<cblas_dgbmv>

=item C<cblas_dtrmv>

=item C<cblas_dtbmv>

=item C<cblas_dtpmv>

=item C<cblas_dtrsv>

=item C<cblas_dtbsv>

=item C<cblas_dtpsv>

=item C<cblas_cgemv>

=item C<cblas_cgbmv>

=item C<cblas_ctrmv>

=item C<cblas_ctbmv>

=item C<cblas_ctpmv>

=item C<cblas_ctrsv>

=item C<cblas_ctbsv>

=item C<cblas_ctpsv>

=item C<cblas_zgemv>

=item C<cblas_zgbmv>

=item C<cblas_ztrmv>

=item C<cblas_ztbmv>

=item C<cblas_ztpmv>

=item C<cblas_ztrsv>

=item C<cblas_ztbsv>

=item C<cblas_ztpsv>

=item C<cblas_ssymv>

=item C<cblas_ssbmv>

=item C<cblas_sspmv>

=item C<cblas_sger>

=item C<cblas_ssyr>

=item C<cblas_sspr>

=item C<cblas_ssyr2>

=item C<cblas_sspr2>

=item C<cblas_dsymv>

=item C<cblas_dsbmv>

=item C<cblas_dspmv>

=item C<cblas_dger>

=item C<cblas_dsyr>

=item C<cblas_dspr>

=item C<cblas_dsyr2>

=item C<cblas_dspr2>

=item C<cblas_chemv>

=item C<cblas_chbmv>

=item C<cblas_chpmv>

=item C<cblas_cgeru>

=item C<cblas_cgerc>

=item C<cblas_cher>

=item C<cblas_chpr>

=item C<cblas_cher2>

=item C<cblas_chpr2>

=item C<cblas_zhemv>

=item C<cblas_zhbmv>

=item C<cblas_zhpmv>

=item C<cblas_zgeru>

=item C<cblas_zgerc>

=item C<cblas_zher>

=item C<cblas_zhpr>

=item C<cblas_zher2>

=item C<cblas_zhpr2>

=item C<cblas_sgemm>

=item C<cblas_ssymm>

=item C<cblas_ssyrk>

=item C<cblas_ssyr2k>

=item C<cblas_strmm>

=item C<cblas_strsm>

=item C<cblas_dgemm>

=item C<cblas_dsymm>

=item C<cblas_dsyrk>

=item C<cblas_dsyr2k>

=item C<cblas_dtrmm>

=item C<cblas_dtrsm>

=item C<cblas_cgemm>

=item C<cblas_csymm>

=item C<cblas_csyrk>

=item C<cblas_csyr2k>

=item C<cblas_ctrmm>

=item C<cblas_ctrsm>

=item C<cblas_zgemm>

=item C<cblas_zsymm>

=item C<cblas_zsyrk>

=item C<cblas_zsyr2k>

=item C<cblas_ztrmm>

=item C<cblas_ztrsm>

=item C<cblas_chemm>

=item C<cblas_cherk>

=item C<cblas_cher2k>

=item C<cblas_zhemm>

=item C<cblas_zherk>

=item C<cblas_zher2k>

=item C<cblas_xerbla>

=back

This module also contains the following constants : 

=over 1

=item C<$CblasRowMajor>
               
=item C<$CblasColMajor> 

=item C<$CblasNoTrans > 

=item C<$CblasTrans > 

=item C<$CblasConjTrans > 

=item C<$CblasUpper> 

=item C<$CblasLower> 

=item C<$CblasNonUnit> 
                             
=item C<$CblasUnit >
 
=item C<$CblasLeft >

=item C<$CblasRight >

=back

For more informations on the functions, we refer you to the GSL offcial documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>




=head1 EXAMPLES

=head1 AUTHORS

Jonathan "Duke" Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2011 Jonathan "Duke" Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

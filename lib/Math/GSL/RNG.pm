# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Don't modify this file, modify the SWIG interface instead.

package Math::GSL::RNG;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::RNGc;
bootstrap Math::GSL::RNG;
package Math::GSL::RNG;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::RNG;

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

package Math::GSL::RNG;

*gsl_rng_types_setup = *Math::GSL::RNGc::gsl_rng_types_setup;
*gsl_rng_alloc = *Math::GSL::RNGc::gsl_rng_alloc;
*gsl_rng_memcpy = *Math::GSL::RNGc::gsl_rng_memcpy;
*gsl_rng_clone = *Math::GSL::RNGc::gsl_rng_clone;
*gsl_rng_free = *Math::GSL::RNGc::gsl_rng_free;
*gsl_rng_set = *Math::GSL::RNGc::gsl_rng_set;
*gsl_rng_max = *Math::GSL::RNGc::gsl_rng_max;
*gsl_rng_min = *Math::GSL::RNGc::gsl_rng_min;
*gsl_rng_name = *Math::GSL::RNGc::gsl_rng_name;
*gsl_rng_fread = *Math::GSL::RNGc::gsl_rng_fread;
*gsl_rng_fwrite = *Math::GSL::RNGc::gsl_rng_fwrite;
*gsl_rng_size = *Math::GSL::RNGc::gsl_rng_size;
*gsl_rng_state = *Math::GSL::RNGc::gsl_rng_state;
*gsl_rng_print_state = *Math::GSL::RNGc::gsl_rng_print_state;
*gsl_rng_env_setup = *Math::GSL::RNGc::gsl_rng_env_setup;
*gsl_rng_get = *Math::GSL::RNGc::gsl_rng_get;
*gsl_rng_uniform = *Math::GSL::RNGc::gsl_rng_uniform;
*gsl_rng_uniform_pos = *Math::GSL::RNGc::gsl_rng_uniform_pos;
*gsl_rng_uniform_int = *Math::GSL::RNGc::gsl_rng_uniform_int;

############# Class : Math::GSL::RNG::gsl_rng_type ##############

package Math::GSL::RNG::gsl_rng_type;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::RNG );
%OWNER = ();
%ITERATORS = ();
*swig_name_get = *Math::GSL::RNGc::gsl_rng_type_name_get;
*swig_name_set = *Math::GSL::RNGc::gsl_rng_type_name_set;
*swig_max_get = *Math::GSL::RNGc::gsl_rng_type_max_get;
*swig_max_set = *Math::GSL::RNGc::gsl_rng_type_max_set;
*swig_min_get = *Math::GSL::RNGc::gsl_rng_type_min_get;
*swig_min_set = *Math::GSL::RNGc::gsl_rng_type_min_set;
*swig_size_get = *Math::GSL::RNGc::gsl_rng_type_size_get;
*swig_size_set = *Math::GSL::RNGc::gsl_rng_type_size_set;
*swig_set_get = *Math::GSL::RNGc::gsl_rng_type_set_get;
*swig_set_set = *Math::GSL::RNGc::gsl_rng_type_set_set;
*swig_get_get = *Math::GSL::RNGc::gsl_rng_type_get_get;
*swig_get_set = *Math::GSL::RNGc::gsl_rng_type_get_set;
*swig_get_double_get = *Math::GSL::RNGc::gsl_rng_type_get_double_get;
*swig_get_double_set = *Math::GSL::RNGc::gsl_rng_type_get_double_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::RNGc::new_gsl_rng_type(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::RNGc::delete_gsl_rng_type($self);
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


############# Class : Math::GSL::RNG::gsl_rng ##############

package Math::GSL::RNG::gsl_rng;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::RNG );
%OWNER = ();
%ITERATORS = ();
*swig_type_get = *Math::GSL::RNGc::gsl_rng_type_get;
*swig_type_set = *Math::GSL::RNGc::gsl_rng_type_set;
*swig_state_get = *Math::GSL::RNGc::gsl_rng_state_get;
*swig_state_set = *Math::GSL::RNGc::gsl_rng_state_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::RNGc::new_gsl_rng(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::RNGc::delete_gsl_rng($self);
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

package Math::GSL::RNG;

*GSL_MAJOR_VERSION = *Math::GSL::RNGc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::RNGc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::RNGc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::RNGc::GSL_NEGZERO;

my %__gsl_rng_borosh13_hash;
tie %__gsl_rng_borosh13_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_borosh13;
$gsl_rng_borosh13= \%__gsl_rng_borosh13_hash;
bless $gsl_rng_borosh13, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_coveyou_hash;
tie %__gsl_rng_coveyou_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_coveyou;
$gsl_rng_coveyou= \%__gsl_rng_coveyou_hash;
bless $gsl_rng_coveyou, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_cmrg_hash;
tie %__gsl_rng_cmrg_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_cmrg;
$gsl_rng_cmrg= \%__gsl_rng_cmrg_hash;
bless $gsl_rng_cmrg, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_fishman18_hash;
tie %__gsl_rng_fishman18_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_fishman18;
$gsl_rng_fishman18= \%__gsl_rng_fishman18_hash;
bless $gsl_rng_fishman18, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_fishman20_hash;
tie %__gsl_rng_fishman20_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_fishman20;
$gsl_rng_fishman20= \%__gsl_rng_fishman20_hash;
bless $gsl_rng_fishman20, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_fishman2x_hash;
tie %__gsl_rng_fishman2x_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_fishman2x;
$gsl_rng_fishman2x= \%__gsl_rng_fishman2x_hash;
bless $gsl_rng_fishman2x, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_gfsr4_hash;
tie %__gsl_rng_gfsr4_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_gfsr4;
$gsl_rng_gfsr4= \%__gsl_rng_gfsr4_hash;
bless $gsl_rng_gfsr4, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_knuthran_hash;
tie %__gsl_rng_knuthran_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_knuthran;
$gsl_rng_knuthran= \%__gsl_rng_knuthran_hash;
bless $gsl_rng_knuthran, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_knuthran2_hash;
tie %__gsl_rng_knuthran2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_knuthran2;
$gsl_rng_knuthran2= \%__gsl_rng_knuthran2_hash;
bless $gsl_rng_knuthran2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_knuthran2002_hash;
tie %__gsl_rng_knuthran2002_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_knuthran2002;
$gsl_rng_knuthran2002= \%__gsl_rng_knuthran2002_hash;
bless $gsl_rng_knuthran2002, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_lecuyer21_hash;
tie %__gsl_rng_lecuyer21_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_lecuyer21;
$gsl_rng_lecuyer21= \%__gsl_rng_lecuyer21_hash;
bless $gsl_rng_lecuyer21, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_minstd_hash;
tie %__gsl_rng_minstd_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_minstd;
$gsl_rng_minstd= \%__gsl_rng_minstd_hash;
bless $gsl_rng_minstd, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_mrg_hash;
tie %__gsl_rng_mrg_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_mrg;
$gsl_rng_mrg= \%__gsl_rng_mrg_hash;
bless $gsl_rng_mrg, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_mt19937_hash;
tie %__gsl_rng_mt19937_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_mt19937;
$gsl_rng_mt19937= \%__gsl_rng_mt19937_hash;
bless $gsl_rng_mt19937, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_mt19937_1999_hash;
tie %__gsl_rng_mt19937_1999_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_mt19937_1999;
$gsl_rng_mt19937_1999= \%__gsl_rng_mt19937_1999_hash;
bless $gsl_rng_mt19937_1999, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_mt19937_1998_hash;
tie %__gsl_rng_mt19937_1998_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_mt19937_1998;
$gsl_rng_mt19937_1998= \%__gsl_rng_mt19937_1998_hash;
bless $gsl_rng_mt19937_1998, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_r250_hash;
tie %__gsl_rng_r250_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_r250;
$gsl_rng_r250= \%__gsl_rng_r250_hash;
bless $gsl_rng_r250, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ran0_hash;
tie %__gsl_rng_ran0_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ran0;
$gsl_rng_ran0= \%__gsl_rng_ran0_hash;
bless $gsl_rng_ran0, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ran1_hash;
tie %__gsl_rng_ran1_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ran1;
$gsl_rng_ran1= \%__gsl_rng_ran1_hash;
bless $gsl_rng_ran1, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ran2_hash;
tie %__gsl_rng_ran2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ran2;
$gsl_rng_ran2= \%__gsl_rng_ran2_hash;
bless $gsl_rng_ran2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ran3_hash;
tie %__gsl_rng_ran3_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ran3;
$gsl_rng_ran3= \%__gsl_rng_ran3_hash;
bless $gsl_rng_ran3, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_rand_hash;
tie %__gsl_rng_rand_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_rand;
$gsl_rng_rand= \%__gsl_rng_rand_hash;
bless $gsl_rng_rand, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_rand48_hash;
tie %__gsl_rng_rand48_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_rand48;
$gsl_rng_rand48= \%__gsl_rng_rand48_hash;
bless $gsl_rng_rand48, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random128_bsd_hash;
tie %__gsl_rng_random128_bsd_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random128_bsd;
$gsl_rng_random128_bsd= \%__gsl_rng_random128_bsd_hash;
bless $gsl_rng_random128_bsd, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random128_glibc2_hash;
tie %__gsl_rng_random128_glibc2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random128_glibc2;
$gsl_rng_random128_glibc2= \%__gsl_rng_random128_glibc2_hash;
bless $gsl_rng_random128_glibc2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random128_libc5_hash;
tie %__gsl_rng_random128_libc5_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random128_libc5;
$gsl_rng_random128_libc5= \%__gsl_rng_random128_libc5_hash;
bless $gsl_rng_random128_libc5, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random256_bsd_hash;
tie %__gsl_rng_random256_bsd_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random256_bsd;
$gsl_rng_random256_bsd= \%__gsl_rng_random256_bsd_hash;
bless $gsl_rng_random256_bsd, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random256_glibc2_hash;
tie %__gsl_rng_random256_glibc2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random256_glibc2;
$gsl_rng_random256_glibc2= \%__gsl_rng_random256_glibc2_hash;
bless $gsl_rng_random256_glibc2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random256_libc5_hash;
tie %__gsl_rng_random256_libc5_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random256_libc5;
$gsl_rng_random256_libc5= \%__gsl_rng_random256_libc5_hash;
bless $gsl_rng_random256_libc5, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random32_bsd_hash;
tie %__gsl_rng_random32_bsd_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random32_bsd;
$gsl_rng_random32_bsd= \%__gsl_rng_random32_bsd_hash;
bless $gsl_rng_random32_bsd, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random32_glibc2_hash;
tie %__gsl_rng_random32_glibc2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random32_glibc2;
$gsl_rng_random32_glibc2= \%__gsl_rng_random32_glibc2_hash;
bless $gsl_rng_random32_glibc2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random32_libc5_hash;
tie %__gsl_rng_random32_libc5_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random32_libc5;
$gsl_rng_random32_libc5= \%__gsl_rng_random32_libc5_hash;
bless $gsl_rng_random32_libc5, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random64_bsd_hash;
tie %__gsl_rng_random64_bsd_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random64_bsd;
$gsl_rng_random64_bsd= \%__gsl_rng_random64_bsd_hash;
bless $gsl_rng_random64_bsd, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random64_glibc2_hash;
tie %__gsl_rng_random64_glibc2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random64_glibc2;
$gsl_rng_random64_glibc2= \%__gsl_rng_random64_glibc2_hash;
bless $gsl_rng_random64_glibc2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random64_libc5_hash;
tie %__gsl_rng_random64_libc5_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random64_libc5;
$gsl_rng_random64_libc5= \%__gsl_rng_random64_libc5_hash;
bless $gsl_rng_random64_libc5, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random8_bsd_hash;
tie %__gsl_rng_random8_bsd_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random8_bsd;
$gsl_rng_random8_bsd= \%__gsl_rng_random8_bsd_hash;
bless $gsl_rng_random8_bsd, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random8_glibc2_hash;
tie %__gsl_rng_random8_glibc2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random8_glibc2;
$gsl_rng_random8_glibc2= \%__gsl_rng_random8_glibc2_hash;
bless $gsl_rng_random8_glibc2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random8_libc5_hash;
tie %__gsl_rng_random8_libc5_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random8_libc5;
$gsl_rng_random8_libc5= \%__gsl_rng_random8_libc5_hash;
bless $gsl_rng_random8_libc5, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random_bsd_hash;
tie %__gsl_rng_random_bsd_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random_bsd;
$gsl_rng_random_bsd= \%__gsl_rng_random_bsd_hash;
bless $gsl_rng_random_bsd, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random_glibc2_hash;
tie %__gsl_rng_random_glibc2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random_glibc2;
$gsl_rng_random_glibc2= \%__gsl_rng_random_glibc2_hash;
bless $gsl_rng_random_glibc2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_random_libc5_hash;
tie %__gsl_rng_random_libc5_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_random_libc5;
$gsl_rng_random_libc5= \%__gsl_rng_random_libc5_hash;
bless $gsl_rng_random_libc5, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_randu_hash;
tie %__gsl_rng_randu_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_randu;
$gsl_rng_randu= \%__gsl_rng_randu_hash;
bless $gsl_rng_randu, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranf_hash;
tie %__gsl_rng_ranf_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranf;
$gsl_rng_ranf= \%__gsl_rng_ranf_hash;
bless $gsl_rng_ranf, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranlux_hash;
tie %__gsl_rng_ranlux_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranlux;
$gsl_rng_ranlux= \%__gsl_rng_ranlux_hash;
bless $gsl_rng_ranlux, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranlux389_hash;
tie %__gsl_rng_ranlux389_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranlux389;
$gsl_rng_ranlux389= \%__gsl_rng_ranlux389_hash;
bless $gsl_rng_ranlux389, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranlxd1_hash;
tie %__gsl_rng_ranlxd1_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranlxd1;
$gsl_rng_ranlxd1= \%__gsl_rng_ranlxd1_hash;
bless $gsl_rng_ranlxd1, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranlxd2_hash;
tie %__gsl_rng_ranlxd2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranlxd2;
$gsl_rng_ranlxd2= \%__gsl_rng_ranlxd2_hash;
bless $gsl_rng_ranlxd2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranlxs0_hash;
tie %__gsl_rng_ranlxs0_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranlxs0;
$gsl_rng_ranlxs0= \%__gsl_rng_ranlxs0_hash;
bless $gsl_rng_ranlxs0, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranlxs1_hash;
tie %__gsl_rng_ranlxs1_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranlxs1;
$gsl_rng_ranlxs1= \%__gsl_rng_ranlxs1_hash;
bless $gsl_rng_ranlxs1, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranlxs2_hash;
tie %__gsl_rng_ranlxs2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranlxs2;
$gsl_rng_ranlxs2= \%__gsl_rng_ranlxs2_hash;
bless $gsl_rng_ranlxs2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_ranmar_hash;
tie %__gsl_rng_ranmar_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_ranmar;
$gsl_rng_ranmar= \%__gsl_rng_ranmar_hash;
bless $gsl_rng_ranmar, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_slatec_hash;
tie %__gsl_rng_slatec_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_slatec;
$gsl_rng_slatec= \%__gsl_rng_slatec_hash;
bless $gsl_rng_slatec, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_taus_hash;
tie %__gsl_rng_taus_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_taus;
$gsl_rng_taus= \%__gsl_rng_taus_hash;
bless $gsl_rng_taus, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_taus2_hash;
tie %__gsl_rng_taus2_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_taus2;
$gsl_rng_taus2= \%__gsl_rng_taus2_hash;
bless $gsl_rng_taus2, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_taus113_hash;
tie %__gsl_rng_taus113_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_taus113;
$gsl_rng_taus113= \%__gsl_rng_taus113_hash;
bless $gsl_rng_taus113, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_transputer_hash;
tie %__gsl_rng_transputer_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_transputer;
$gsl_rng_transputer= \%__gsl_rng_transputer_hash;
bless $gsl_rng_transputer, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_tt800_hash;
tie %__gsl_rng_tt800_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_tt800;
$gsl_rng_tt800= \%__gsl_rng_tt800_hash;
bless $gsl_rng_tt800, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_uni_hash;
tie %__gsl_rng_uni_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_uni;
$gsl_rng_uni= \%__gsl_rng_uni_hash;
bless $gsl_rng_uni, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_uni32_hash;
tie %__gsl_rng_uni32_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_uni32;
$gsl_rng_uni32= \%__gsl_rng_uni32_hash;
bless $gsl_rng_uni32, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_vax_hash;
tie %__gsl_rng_vax_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_vax;
$gsl_rng_vax= \%__gsl_rng_vax_hash;
bless $gsl_rng_vax, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_waterman14_hash;
tie %__gsl_rng_waterman14_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_waterman14;
$gsl_rng_waterman14= \%__gsl_rng_waterman14_hash;
bless $gsl_rng_waterman14, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_zuf_hash;
tie %__gsl_rng_zuf_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_zuf;
$gsl_rng_zuf= \%__gsl_rng_zuf_hash;
bless $gsl_rng_zuf, Math::GSL::RNG::gsl_rng_type;

my %__gsl_rng_default_hash;
tie %__gsl_rng_default_hash,"Math::GSL::RNG::gsl_rng_type", $Math::GSL::RNGc::gsl_rng_default;
$gsl_rng_default= \%__gsl_rng_default_hash;
bless $gsl_rng_default, Math::GSL::RNG::gsl_rng_type;
*gsl_rng_default_seed = *Math::GSL::RNGc::gsl_rng_default_seed;

@EXPORT_OK = qw/
                 gsl_rng_alloc gsl_rng_set gsl_rng_get gsl_rng_free gsl_rng_memcpy
                 gsl_rng_fwrite gsl_rng_fread gsl_rng_clone gsl_rng_max gsl_rng_min
                 gsl_rng_name gsl_rng_size gsl_rng_state gsl_rng_print_state gsl_rng_uniform gsl_rng_uniform_pos gsl_rng_uniform_int 
                $gsl_rng_default $gsl_rng_knuthran $gsl_rng_ran0 $gsl_rng_borosh13
                $gsl_rng_coveyou $gsl_rng_cmrg $gsl_rng_fishman18 $gsl_rng_fishman20 $gsl_rng_fishman2x
                $gsl_rng_gfsr4 $gsl_rng_knuthran $gsl_rng_knuthran2 $gsl_rng_knuthran2002 $gsl_rng_lecuyer21
                $gsl_rng_minstd $gsl_rng_mrg $gsl_rng_mt19937 $gsl_rng_mt19937_1999 $gsl_rng_mt19937_1998
                $gsl_rng_r250 $gsl_rng_ran0 $gsl_rng_ran1 $gsl_rng_ran2 $gsl_rng_ran3
                $gsl_rng_rand $gsl_rng_rand48 $gsl_rng_random128_bsd $gsl_rng_random128_gli $gsl_rng_random128_lib
                $gsl_rng_random256_bsd $gsl_rng_random256_gli $gsl_rng_random256_lib $gsl_rng_random32_bsd
                $gsl_rng_random32_glib $gsl_rng_random32_libc $gsl_rng_random64_bsd $gsl_rng_random64_glib 
                $gsl_rng_random64_libc $gsl_rng_random8_bsd $gsl_rng_random8_glibc $gsl_rng_random8_libc5
                $gsl_rng_random_bsd $gsl_rng_random_glibc2 $gsl_rng_random_libc5 $gsl_rng_randu
                $gsl_rng_ranf $gsl_rng_ranlux $gsl_rng_ranlux389 $gsl_rng_ranlxd1 $gsl_rng_ranlxd2 $gsl_rng_ranlxs0
                $gsl_rng_ranlxs1 $gsl_rng_ranlxs2 $gsl_rng_ranmar $gsl_rng_slatec $gsl_rng_taus $gsl_rng_taus2
                $gsl_rng_taus113 $gsl_rng_transputer $gsl_rng_tt800 $gsl_rng_uni $gsl_rng_uni32 $gsl_rng_vax
                $gsl_rng_waterman14 $gsl_rng_zuf 
              /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

=head1 NAME

Math::GSL::RNG - Random Number Generators

=head1 SYNOPSIS

    use Math::GSL::RNG;
    my $rng     = Math::GSL::RNG->new;
    my @random  = $rng->get(100);

=head2 Math::GSL::RNG->new($type, $seed)

    my $rng = Math::GSL::RNG->new;
    my $rng = Math::GSL::RNG->new($gsl_rng_knuthran,5);

Creates a new RNG object of type $type, seeded with $seed. Both of these
parameters are optional. The type $gsl_rng_default is used when no $type
is given.

=cut
 
sub new {
    my ($class, $type, $seed) = @_;
    $type ||= $gsl_rng_default;
    $seed ||= int 100*rand;

    my $self = {};
    my $rng  = gsl_rng_alloc($type);
    gsl_rng_set($rng, $seed);

    $self->{_rng} = $rng; 
    bless $self, $class;
}

=head2 copy()

    my $copy = $rng->copy;

Make a copy of a RNG object.

=cut

sub copy {
    my ($self)    = @_;
    my $copy      = Math::GSL::RNG->new;
    $copy->{_rng} = gsl_rng_clone($self->{_rng});

    return $copy;
}

=head2 free()

    $rng->free();

Free memory associated with RNG object.

=cut

sub free {
    my ($self)    = @_;
    gsl_rng_free($self->{_rng});
}

=head2 name()

   my $name = $rng->name();

Get the name of the RNG object as a string.

=cut

sub name {
    my ($self)    = @_;
    gsl_rng_name($self->{_rng});
}

=head2 get()

    my $nextval  = $rng->get;
    my (@values) = $rng->get(100);

Get the next random value from the RNG object. If given an integer N, returns the next N values.

=cut

sub get {
    my ($self, $num_values) = @_;
    $num_values ||= 1;

    return map { gsl_rng_get($self->{_rng}) } (1 .. $num_values);
}

=head2 raw()

    my $raw = $rng->raw();

Return the raw GSL RNG object, useful for functions which take a RNG, such as the Monte Carlo integration functions or the random number distribution functions in Math::GSL::Randist.

=cut

sub raw {
    my $self = shift;
    return $self->{_rng};
}

__END__


=head1 DESCRIPTION

=over 1

=item gsl_rng_alloc($T) - This function returns a pointer to a newly-created instance of a random number generator of type $T. $T must be one of the constants below. The generator is automatically initialized with the default seed, $gsl_rng_default.

=item gsl_rng_set($r, $s) - This function initializes (or `seeds') the random number generator. If the generator is seeded with the same value of $s on two different runs, the same stream of random numbers will be generated by successive calls to the routines below. If different values of $s are supplied, then the generated streams of random numbers should be completely different. If the seed $s is zero then the standard seed from the original implementation is used instead. For example, the original Fortran source code for the ranlux generator used a seed of 314159265, and so choosing $s equal to zero reproduces this when using $gsl_rng_ranlux.

=item gsl_rng_get($r) - This function returns a random integer from the generator $r. The minimum and maximum values depend on the algorithm used, but all integers in the range [min,max] are equally likely. The values of min and max can determined using the auxiliary functions gsl_rng_max($r) and gsl_rng_min($r).

=item gsl_rng_free($r) - This function frees all the memory associated with the generator $r.

=item gsl_rng_memcpy($dest, $src) - This function copies the random number generator $src into the pre-existing generator $dest, making $dest into an exact copy of $src. The two generators must be of the same type.

=item gsl_rng_uniform($r) - This function returns a double precision floating point number uniformly distributed in the range [0,1). The range includes 0.0 but excludes 1.0. The value is typically obtained by dividing the result of gsl_rng_get($r) by gsl_rng_max($r) + 1.0 in double precision. Some generators compute this ratio internally so that they can provide floating point numbers with more than 32 bits of randomness (the maximum number of bits that can be portably represented in a single unsigned long int). 

=item gsl_rng_uniform_pos($r) - This function returns a positive double precision floating point number uniformly distributed in the range (0,1), excluding both 0.0 and 1.0. The number is obtained by sampling the generator with the algorithm of gsl_rng_uniform until a non-zero value is obtained. You can use this function if you need to avoid a singularity at 0.0.

=item gsl_rng_uniform_int($r, $n) - This function returns a random integer from 0 to $n-1 inclusive by scaling down and/or discarding samples from the generator $r. All integers in the range [0,$n-1] are produced with equal probability. For generators with a non-zero minimum value an offset is applied so that zero is returned with the correct probability. Note that this function is designed for sampling from ranges smaller than the range of the underlying generator. The parameter $n must be less than or equal to the range of the generator $r. If $n is larger than the range of the generator then the function calls the error handler with an error code of $GSL_EINVAL and returns zero. In particular, this function is not intended for generating the full range of unsigned integer values [0,2^32-1]. Instead choose a generator with the maximal integer range and zero mimimum value, such as $gsl_rng_ranlxd1, $gsl_rng_mt19937 or $gsl_rng_taus, and sample it directly using gsl_rng_get. The range of each generator can be found using the auxiliary functions described in the next section. 

=item gsl_rng_fwrite($stream, $r) - This function writes the random number state of the random number generator $r to the stream $stream (opened with the gsl_fopen function from the Math::GSL module) in binary format. The return value is 0 for success and $GSL_EFAILED if there was a problem writing to the file. Since the data is written in the native binary format it may not be portable between different architectures.

=item gsl_rng_fread($stream, $r) - This function reads the random number state into the random number generator $r from the open stream $stream (opened with the gsl_fopen function from the Math::GSL module) in binary format. The random number generator $r must be preinitialized with the correct random number generator type since type information is not saved. The return value is 0 for success and $GSL_EFAILED if there was a problem reading from the file. The data is assumed to have been written in the native binary format on the same architecture.

=item gsl_rng_clone($r) - This function returns a pointer to a newly created generator which is an exact copy of the generator $r.

=item gsl_rng_max($r) - This function returns the largest value that gsl_rng_get can return.

=item gsl_rng_min($r) - gsl_rng_min returns the smallest value that gsl_rng_get can return. Usually this value is zero. There are some generators with algorithms that cannot return zero, and for these generators the minimum value is 1.

=item gsl_rng_name($r) - This function returns a pointer to the name of the generator. For example,

=over

=item print "r is a " . gsl_rng_name($r) . "generator\n";

=item would print something like r is a 'taus' generator. 

=back

=item gsl_rng_size($r) - This function returns the size of the state of generator $r. You can use this information to access the state directly.

=item gsl_rng_state($r) - This function returns a pointer to the state of generator $r. You can use this information to access the state directly.

=item gsl_rng_print_state($r) 

=back

=head1 Random Number Generator Types

=over 1

=item $gsl_rng_default

=item $gsl_rng_knuthran

=item $gsl_rng_ran0

=item $gsl_rng_borosh13

=item $gsl_rng_coveyou

=item $gsl_rng_cmrg

=item $gsl_rng_fishman18

=item $gsl_rng_fishman20

=item $gsl_rng_fishman2x - This is the L'Ecuyer–Fishman random number generator. It is taken from Knuth's Seminumerical Algorithms, 3rd Ed., page 108. Its sequence is, z_{n+1} = (x_n - y_n) mod m with m = 2^31 - 1. x_n and y_n are given by the fishman20 and lecuyer21 algorithms. The seed specifies the initial value, x_1. 

=item $gsl_rng_gfsr4

=item $gsl_rng_knuthran

=item $gsl_rng_knuthran2

=item $gsl_rng_knuthran2002

=item $gsl_rng_lecuyer21

=item $gsl_rng_minstd

=item $gsl_rng_mrg

=item $gsl_rng_mt19937

=item $gsl_rng_mt19937_1999

=item $gsl_rng_mt19937_1998

=item $gsl_rng_r250

=item $gsl_rng_ran0

=item $gsl_rng_ran1

=item $gsl_rng_ran2

=item $gsl_rng_ran3

=item $gsl_rng_rand - This is the BSD rand generator. Its sequence is x_{n+1} = (a x_n + c) mod m with a = 1103515245, c = 12345 and m = 2^31. The seed specifies the initial value, x_1. The period of this generator is 2^31, and it uses 1 word of storage per generator. 

=item $gsl_rng_rand48

=item $gsl_rng_random128_bsd

=item $gsl_rng_random128_gli

=item $gsl_rng_random128_lib

=item $gsl_rng_random256_bsd

=item $gsl_rng_random256_gli

=item $gsl_rng_random256_lib

=item $gsl_rng_random32_bsd

=item $gsl_rng_random32_glib

=item $gsl_rng_random32_libc

=item $gsl_rng_random64_bsd

=item $gsl_rng_random64_glib

=item $gsl_rng_random64_libc

=item $gsl_rng_random8_bsd

=item $gsl_rng_random8_glibc

=item $gsl_rng_random8_libc5

=item $gsl_rng_random_bsd

=item $gsl_rng_random_glibc2

=item $gsl_rng_random_libc5

=item $gsl_rng_randu

=item $gsl_rng_ranf

=item $gsl_rng_ranlux

=item $gsl_rng_ranlux389

=item $gsl_rng_ranlxd1

=item $gsl_rng_ranlxd2

=item $gsl_rng_ranlxs0

=item $gsl_rng_ranlxs1

=item $gsl_rng_ranlxs2

=item $gsl_rng_ranmar - This is the RANMAR lagged-fibonacci generator of Marsaglia, Zaman and Tsang. It is a 24-bit generator, originally designed for single-precision IEEE floating point numbers. It was included in the CERNLIB high-energy physics library.

=item $gsl_rng_slatec - This is the SLATEC random number generator RAND. It is ancient. The original source code is available from NETLIB.

=item $gsl_rng_taus

=item $gsl_rng_taus2

=item $gsl_rng_taus113

=item $gsl_rng_transputer

=item $gsl_rng_tt800

=item $gsl_rng_uni

=item $gsl_rng_uni32

=item $gsl_rng_vax - This is the VAX generator MTH$RANDOM. Its sequence is, x_{n+1} = (a x_n + c) mod m with a = 69069, c = 1 and m = 2^32. The seed specifies the initial value, x_1. The period of this generator is 2^32 and it uses 1 word of storage per generator. 

=item $gsl_rng_waterman14

=item $gsl_rng_zuf - This is the ZUFALL lagged Fibonacci series generator of Peterson. Its sequence is,

=over

=item          t = u_{n-273} + u_{n-607}

=item          u_n  = t - floor(t)
     
=back

 The original source code is available from NETLIB. For more information see,

 * W. Petersen, “Lagged Fibonacci Random Number Generators for the NEC SX-3”, International Journal of High Speed Computing (1994). 

=back

For more informations on the functions, we refer you to the GSL offcial documentation: 

L<http://www.gnu.org/software/gsl/manual/html_node/>

Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want

=head1 EXAMPLES

The following example will print out a list a random integers between certain
minimum and maximum values. The command line arguments are first the number of
random numbers wanted, the minimum and then maximum. The defaults are 10, 0 and
100, respectively. 

    use Math::GSL::RNG qw/:all/;
    my $seed = int rand(100);
    my $rng  = Math::GSL::RNG->new($gsl_rng_knuthran, $seed );
    my ($num,$min,$max) = @ARGV;
    $num ||= 10;
    $min ||= 0;
    $max ||= 100;
    print join "\n", map { $min + $rng->get % ($max-$min+1)  } (1..$num);
    print "\n";

The C<$seed> argument is optional but encouraged. This program is available in
the B<examples/> directory that comes with the source of this module.

If you would like a series of random non-integer numbers, then you can generate one "scaling factor" 
and multiple by that, such as

    use Math::GSL::RNG qw/:all/;
    my $scale= rand(10);
    my $seed = int rand(100);
    my $rng  = Math::GSL::RNG->new($gsl_rng_knuthran, $seed );
    my ($num,$min,$max) = (10,0,100);
    print join "\n", map { $scale*($min + $rng->get % ($max-$min+1))  } (1..$num);
    print "\n";

=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;

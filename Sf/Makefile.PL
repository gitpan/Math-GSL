use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

my $gsl_inc = `gsl-config --cflags`;
chomp $gsl_inc;

my $gsl_libs = `gsl-config --libs`;
chomp $gsl_libs;


WriteMakefile(
    'NAME'      => 'Math::Gsl::Sf',
    'VERSION_FROM' => 'Sf.pm', # finds $VERSION
    'LIBS'              => [ $gsl_libs ], 
    'DEFINE'            => '', 
    'INC'               => $gsl_inc ,

);


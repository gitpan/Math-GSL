use ExtUtils::MakeMaker;

my $gsl_ver = `gsl-config --version`;
chomp $gsl_ver;
die "Sorry, Gsl >= v0.94 is required" if $gsl_ver < 0.94;

my $gsl_inc = `gsl-config --cflags`;
chomp $gsl_inc;

my $gsl_libs = `gsl-config --libs`;
chomp $gsl_libs;

WriteMakefile(
    'DIR'  		=> [qw(Sf Polynomial)],
    'NAME'		=> 'Math::Gsl',
    'VERSION_FROM'	=> 'Gsl.pm', 
    'PREREQ_PM'		=> {}, 
    ($] >= 5.005 ?    
      (ABSTRACT => 'Perl Interface to GSL',
       AUTHOR     => 'Jonathan Leto <jonathan@leto.net>') : ()),
    'LIBS'		=> [$gsl_libs], 
    'DEFINE'		=> '', 
    'INC'		=> $gsl_inc,
);


use ExtUtils::MakeMaker;

my $gsl_ver = `gsl-config --version`;
chomp $gsl_ver;
die "Sorry, Gsl >= v0.94 is required" if $gsl_ver < 0.94;

WriteMakefile(
    'DIR'  		=> [qw(Sf)],
    'NAME'		=> 'Math::Gsl',
    'VERSION_FROM'	=> 'Gsl.pm', 
    'PREREQ_PM'		=> {}, 
    ($] >= 5.005 ?    
      (ABSTRACT => 'Perl Interface to GSL',
       AUTHOR     => 'Jonathan Leto <jonathan@leto.net>') : ()),
    'LIBS'		=> ['-lgsl -lgslcblas -lm'], 
    'DEFINE'		=> '', 
    'INC'		=> '', 
);


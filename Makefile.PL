use ExtUtils::MakeMaker;

WriteMakefile(
  NAME => 'Math::GSL',
  VERSION_FROM => 'GSL.pm',
  LIBS => ['-lgsl -lgslcblas -lm'],
  );

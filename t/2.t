use Test::More (tests => 21);
use Math::Complex;
use Math::Gsl::Polynomial qw(poly_complex_solve);
ok(1);

# BEGIN { use_ok('Math::Gsl::Polynomial qw(poly_complex_solve)') }

my @a = (15,-8,1);

my @roots = poly_complex_solve(@a);
print "#@roots\n";
is(4,scalar(@roots),"Solved");
is($roots[0],3,"1st root is 3");
is($roots[1],0,"1st root is real");
is($roots[2],5,"2nd root is 5");
is($roots[3],0,"2nd root is real");

@roots = poly_complex_solve(-1,0,0,0,0,1);
print "#@roots\n";
is(scalar(@roots),10,"Solved");


while (@roots)
 {
  my $n = Math::Complex->new(splice(@roots,0,2));
  ok(abs($n**5-1) < 1.0e-6,"Fifth root");
 }

my @poly = (1,-3.17771244049072,3.9795618057251,-1.72559440135956,-0.857469737529755,0.766406536102295,0.816600441932678,-1.16691339015961,0.416294276714325);

@roots = poly_complex_solve(@poly);
print "#@roots\n";
is(scalar(@roots),16,"Solved");

while (@roots)
 {
  my $n = Math::Complex->new(splice(@roots,0,2));
  my $v = poly_complex_eval($n,@poly);
  print "# $n => $v\n";
  ok(abs($v) < 1.0e-6,"Is a solution");
 }


sub poly_complex_eval
{
 my ($n,@a) = @_;
 my $v = 0;
 while (@a)
  {
   $v = $v*$n + pop(@a);
  }
 return $v;
}

	     
	     
	     
	     

use Test;
BEGIN{  plan tests => 10; }

use strict;
my $eps = .0000000000001;

use Math::Gsl;
ok(1);
use Math::Gsl::Sf qw(gamma);
ok(1);

#########################

#check if we can use fully qualified functions
ok( equals(Math::Gsl::Sf::gamma(5), 24) );

#check if we can use direct XS func ( no wrapper, much quicker )
ok( equals(Math::Gsl::Sf::sf_gamma(5), 24) );

#check if imported func is ok
ok( equals(gamma(5), 24) );

#check if we can create a new object, and it is a hashref
my $sf = new Math::Gsl::Sf;
ok(  $sf->isa("Math::Gsl::Sf")  );

#check if we can use the object
ok( equals($sf->gamma(5), 24) );

#check if we can create a new result obj, and it is a scalar ref
my $r = new Math::Gsl::Sf::Result;
ok( ref $r eq "gsl_sf_resultPtr" );

#check with using result obj
my $status = $sf->gamma_e( 5 , $r );
ok( equals($r->val , 24) );

#check fully qualified with result obj
$status = Math::Gsl::Sf::gamma_e(5, $r);
ok( equals($r->val, 24) );

### das boot

# take 2 values, and see if they
# are equal to a certain tolerance
sub equals {
	my ($x,$y) = @_;
	return 1 if ( abs( $x - $y) < $eps );
	return 0;
}
	
	


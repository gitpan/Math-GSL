# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl 1.t'

#########################

use Test;
use strict;

BEGIN { plan tests => 10 };
use Math::Gsl;
ok(1); 
use Math::Gsl::Sf qw(gamma);
ok(1);

#########################

#check if we can use fully qualified functions
ok( Math::Gsl::Sf::gamma(5), 24);

#check if we can use direct XS func ( no wrapper, much quicker )
ok( Math::Gsl::Sf::sf_gamma(5), 24 );

#check if imported func is ok
ok( gamma(5), 24 );

#check if we can create a new object, and it is a hashref
my $sf = new Math::Gsl::Sf;
ok( ref $sf, "Math::Gsl::Sf" ); 

#check if we can use the object
ok( $sf->gamma(5), 24 );

#check if we can create a new result obj, and it is a scalar ref
my $r = new Math::Gsl::Sf::Result;
ok( ref $r, "gsl_sf_resultPtr" );

#check with using result obj
my $status = $sf->gamma_e( 5 , $r );
ok( $r->val , 24 );

#check fully qualified with result obj
$status = Math::Gsl::Sf::gamma_e(5, $r);
ok( $r->val, 24 );




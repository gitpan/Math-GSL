package Math::Gsl::Polynomial;
use Exporter;
use base 'Exporter';
use XSLoader;
use vars qw($VERSION @EXPORT);
$VERSION = '0.08';
@EXPORT = ();
@EXPORT_OK = qw(poly_complex_solve);

XSLoader::load 'Math::Gsl::Polynomial';

1;
__END__

=head1 NAME

Math::Gsl::Polynomial - Perl Interface to solving polynomials with The GNU Scientific Library

=head1 SYNOPSIS

	use Math::Complex;
	use Math::Gsl::Polynomial qw(poly_complex_solve);
	my @a = (15,-8,1);
	my @roots = poly_complex_solve(@a);

=head1 DESCRIPTION

  The only function defined in this module is poly_complex_solve. Note
  that it is not imported by default.

=head1 Exported Constants

  None

=head1 Exported Functions

  poly_complex_solve

=head1 SEE ALSO

GNU Scientific Library http://sources.redhat.com/gsl

To get Math::Gsl: http://www.leto.net/code/gsl/

Math::Polynomial

Math::Polynomial::Solve

=cut

package Math::GSL;
use Exporter;
use base 'Exporter';
use XSLoader;
use vars qw($VERSION @EXPORT);
$VERSION = '0.01';
@EXPORT = qw(poly_complex_solve);

XSLoader::load 'Math::GSL';

1;
__END__

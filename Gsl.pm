package Math::Gsl;

use strict;
no strict 'refs';
use Carp;
use vars qw($VERSION %EXPORT_TAGS @ISA @EXPORT @EXPORT_OK $AUTOLOAD);

require Exporter;
#require AutoLoader;
require DynaLoader;


@ISA = qw(Exporter DynaLoader);
@EXPORT = ();
%EXPORT_TAGS = ();

# maximum smarts, minimum effort
# fills in EXPORT_OK with all functions, then fills in :all tag with
# them all too, so they can all be slurped in at once
for ( keys %EXPORT_TAGS ){
	Exporter::export_ok_tags($_);
	push(@{ $EXPORT_TAGS{"all"} },  @{$EXPORT_TAGS{$_}} );
}

$VERSION = '0.08';

sub AUTOLOAD {
    # This AUTOLOAD is used to 'autoload' constants from the constant()
    # XS function.  If a constant is not found then control is passed
    # to the AUTOLOAD in AutoLoader.

    my $constname;
    ($constname = $AUTOLOAD) =~ s/.*:://;
    croak "& not defined" if $constname eq 'constant';
    my $val = constant($constname, @_ ? $_[0] : 0);
    if ($! != 0) {
	if ($! =~ /Invalid/) {
	    $AutoLoader::AUTOLOAD = $AUTOLOAD;
	    goto &AutoLoader::AUTOLOAD;
	}
	else {
		croak "Your vendor has not defined Math::Gsl macro $constname";
	}
    }
    no strict 'refs';
    *$AUTOLOAD = sub () { $val };
    goto &$AUTOLOAD;
}

bootstrap Math::Gsl $VERSION;


1;
__END__

=head1 NAME

Gsl - Perl Interface to The GNU Scientific Library

=head1 SYNOPSIS

  use Math::Gsl;
  use Math::Gsl::Sf;	     # use Special Function library
  use Math::Gsl::Polynomial; # use poly_complex_solve


=head1 DESCRIPTION

  Currently this module implements the GSL Special function library and the
  single GSL function poly_complex_solve.
  Please see the manual page for Math::Gsl::Sf for relevant documenation.


=head1 Exported constants

	None

=head1 Exported functions

	None

=head1 AUTHOR

Jonathan Leto, jonathan@leto.net

=head1 SEE ALSO

GNU Scientific Library http://sources.redhat.com/gsl

To get Math::Gsl: http://www.leto.net/code/gsl/

=cut

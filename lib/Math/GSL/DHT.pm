# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

package Math::GSL::DHT;
use base qw(Exporter);
use base qw(DynaLoader);
package Math::GSL::DHTc;
bootstrap Math::GSL::DHT;
package Math::GSL::DHT;
@EXPORT = qw();

# ---------- BASE METHODS -------------

package Math::GSL::DHT;

sub TIEHASH {
    my ($classname,$obj) = @_;
    return bless $obj, $classname;
}

sub CLEAR { }

sub FIRSTKEY { }

sub NEXTKEY { }

sub FETCH {
    my ($self,$field) = @_;
    my $member_func = "swig_${field}_get";
    $self->$member_func();
}

sub STORE {
    my ($self,$field,$newval) = @_;
    my $member_func = "swig_${field}_set";
    $self->$member_func($newval);
}

sub this {
    my $ptr = shift;
    return tied(%$ptr);
}


# ------- FUNCTION WRAPPERS --------

package Math::GSL::DHT;

*gsl_dht_alloc = *Math::GSL::DHTc::gsl_dht_alloc;
*gsl_dht_new = *Math::GSL::DHTc::gsl_dht_new;
*gsl_dht_init = *Math::GSL::DHTc::gsl_dht_init;
*gsl_dht_x_sample = *Math::GSL::DHTc::gsl_dht_x_sample;
*gsl_dht_k_sample = *Math::GSL::DHTc::gsl_dht_k_sample;
*gsl_dht_free = *Math::GSL::DHTc::gsl_dht_free;
*gsl_dht_apply = *Math::GSL::DHTc::gsl_dht_apply;

############# Class : Math::GSL::DHT::gsl_dht_struct ##############

package Math::GSL::DHT::gsl_dht_struct;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Math::GSL::DHT );
%OWNER = ();
%ITERATORS = ();
*swig_size_get = *Math::GSL::DHTc::gsl_dht_struct_size_get;
*swig_size_set = *Math::GSL::DHTc::gsl_dht_struct_size_set;
*swig_nu_get = *Math::GSL::DHTc::gsl_dht_struct_nu_get;
*swig_nu_set = *Math::GSL::DHTc::gsl_dht_struct_nu_set;
*swig_xmax_get = *Math::GSL::DHTc::gsl_dht_struct_xmax_get;
*swig_xmax_set = *Math::GSL::DHTc::gsl_dht_struct_xmax_set;
*swig_kmax_get = *Math::GSL::DHTc::gsl_dht_struct_kmax_get;
*swig_kmax_set = *Math::GSL::DHTc::gsl_dht_struct_kmax_set;
*swig_j_get = *Math::GSL::DHTc::gsl_dht_struct_j_get;
*swig_j_set = *Math::GSL::DHTc::gsl_dht_struct_j_set;
*swig_Jjj_get = *Math::GSL::DHTc::gsl_dht_struct_Jjj_get;
*swig_Jjj_set = *Math::GSL::DHTc::gsl_dht_struct_Jjj_set;
*swig_J2_get = *Math::GSL::DHTc::gsl_dht_struct_J2_get;
*swig_J2_set = *Math::GSL::DHTc::gsl_dht_struct_J2_set;
sub new {
    my $pkg = shift;
    my $self = Math::GSL::DHTc::new_gsl_dht_struct(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Math::GSL::DHTc::delete_gsl_dht_struct($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


# ------- VARIABLE STUBS --------

package Math::GSL::DHT;

*GSL_MAJOR_VERSION = *Math::GSL::DHTc::GSL_MAJOR_VERSION;
*GSL_MINOR_VERSION = *Math::GSL::DHTc::GSL_MINOR_VERSION;
*GSL_POSZERO = *Math::GSL::DHTc::GSL_POSZERO;
*GSL_NEGZERO = *Math::GSL::DHTc::GSL_NEGZERO;

@EXPORT_OK = qw/
               gsl_dht_alloc 
               gsl_dht_new 
               gsl_dht_init 
               gsl_dht_x_sample 
               gsl_dht_k_sample 
               gsl_dht_free 
               gsl_dht_apply 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::DHT - Discrete Hankel Transforms

=head1 SYNOPSIS

    use Math::GSL::DHT qw/:all/;

=head1 DESCRIPTION

Here is a list of all the functions included in this module :

=over

=item C<gsl_dht_alloc($size)> - This function allocates a Discrete Hankel transform object of size $size.  

=item C<gsl_dht_new($size, $nu, $xmax)> -  This function allocates a Discrete Hankel transform object of size $size and initializes it for the given values of $nu and $xmax. 

=item C<gsl_dht_init($t, $nu, $xmax)> - This function initializes the transform $t for the given values of $nu and $xmax.  

=item C<gsl_dht_x_sample($t, $n)> - This function returns the value of the $n-th sample point in the unit interval, (j_{\nu,n+1}/j_{\nu,M}) X. These are the points where the function f(t) is assumed to be sampled.  

=item C<gsl_dht_k_sample($t, $n)> - This function returns the value of the $n-th sample point in "k-space", j_{\nu,n+1}/X.  

=item C<gsl_dht_free($t)> - This function frees the transform $t.  

=item C<gsl_dht_apply> 

=back

=head1 EXAMPLES

=head1 AUTHORS

Jonathan "Duke" Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2011 Jonathan "Duke" Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
1;

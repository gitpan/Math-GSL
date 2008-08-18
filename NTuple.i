%module "Math::GSL::NTuple"

%typemap(in) void *ntuple_data {
    fprintf(stderr,"symname=$symname\n");
    if ($input) 
        $1 = (double *) $input;
};

%typemap(argout) void *ntuple_data {
    //Perl_sv_dump($1);
}

%{
    #include "gsl/gsl_ntuple.h"
%}

%include "gsl/gsl_ntuple.h"


%perlcode %{

# Intermittent failure happens *after* this
# END { warn "This is the end" }


@EXPORT_OK = qw/
               gsl_ntuple_open 
               gsl_ntuple_create 
               gsl_ntuple_write 
               gsl_ntuple_read 
               gsl_ntuple_bookdata 
               gsl_ntuple_project 
               gsl_ntuple_close 
             /;
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

__END__

=head1 NAME

Math::GSL::NTuple - Functions for creating and manipulating ntuples, sets of values associated with events

=head1 SYNOPSIS

This module is partially implemented. Patches Welcome!

    use Math::GSL::NTuple qw /:all/;

=head1 DESCRIPTION

Here is a list of all the functions in this module :

=over

=item * <gsl_ntuple_open($filename, $ntuple_data, $size)> - This function opens an existing ntuple file $filename for reading and returns a pointer to a corresponding ntuple struct. The ntuples in the file must have size $size. A pointer to memory for the current ntuple row $ntuple_data, which is an array reference, must be supplied—this is used to copy ntuples in and out of the file.

=item * <gsl_ntuple_create > - This function creates a new write-only ntuple file $filename for ntuples of size $size and returns a pointer to the newly created ntuple struct. Any existing file with the same name is truncated to zero length and overwritten. A pointer to memory for the current ntuple row $ntuple_data, which is an array reference, must be supplied—this is used to copy ntuples in and out of the file. 


=item * <gsl_ntuple_write($ntuple)> - This function writes the current $ntuple $ntuple->{ntuple_data} of size $ntuple->{size} to the corresponding file.

=item * <gsl_ntuple_bookdata($ntuple)> - This function is a synonym for gsl_ntuple_write.

=item * <gsl_ntuple_read($ntuple)> - This function reads the current row of the ntuple file for ntuple and stores the values in $ntuple->{data}.

=item * <gsl_ntuple_project >

=item * <gsl_ntuple_close($ntuple)> - This function closes the ntuple file ntuple and frees its associated allocated memory. 

=back

For more informations on the functions, we refer you to the GSL offcial
documentation: L<http://www.gnu.org/software/gsl/manual/html_node/>

 Tip : search on google: site:http://www.gnu.org/software/gsl/manual/html_node/ name_of_the_function_you_want


=head1 AUTHORS

Jonathan Leto <jonathan@leto.net> and Thierry Moisan <thierry.moisan@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 Jonathan Leto and Thierry Moisan

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

%}

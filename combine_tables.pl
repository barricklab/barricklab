#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

combine_tables.pl

=head1 SYNOPSIS

Usage: combine_tables.pl [-c -f CSV] -o output_file input_file_1 input_file_2 [input_file_3...]

Concatenates text files, leaving off the first line (header) of each file except for the first. 
Checks that the header is identical between all files. Optionally, adds a new column with the
base filename (no extension) that originally contained that data. When adding this new column
assumes CSV format unless TSV is specified at the command line.

Type "perldoc combine_tables.pl" or "combine_tables.pl --man" for more detailed help!

=head1 DESCRIPTION

=over

=item B<-o,--output> <path> 

Only execute for directories/files that match this regular expression.
Names beginning with a period '.' or underscore '_' are always ignored.

=item B<-f,--format> <format> 

Valid choices are CSV and TSV. Only needed when using -c option. (Default: CSV)

=item B<-c,--file-name-column> 

Add a new column with the file name that each line was from to the output. 
All extensions are stripped from the file name.

=back

=head1 AUTHOR

Jeffrey Barrick

=head1 COPYRIGHT

Copyright 2018.  All rights reserved.

=cut

###
# End Pod Documentation
###

use strict;

use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;

#Get options
use Getopt::Long;
use Pod::Usage;
my ($help, $man);
my ($output, $format, $file_name_column);
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'output|o=s' => \$output,
	'format|f=s' => \$format,
	'file-name-column|c' => \$file_name_column,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

$format='CSV' if (!$format);

my @input_file_names = @ARGV;

scalar @input_file_names > 1 or die "Must be used to combine at least two input files";
defined $output or die "Must provide output file name (-o)";
open OUT, ">$output" or die "Could not open output file $output";

my $output_header_line;

foreach my $input_file_name (@input_file_names) {
  open IN, "<$input_file_name" or die "Could not open input file $input_file_name";
  
  my $base_file_name = $input_file_name;
  $base_file_name =~ s/^.+\///;
  $base_file_name =~ s/\..+$//;
  
  my $input_header_line = <IN>;
  if (!defined $output_header_line) {
    
    my $printed_input_header_line = $input_header_line;
    if ($file_name_column) {
         $printed_input_header_line = "file" . ( ($format eq 'CSV') ? "," : "\t") . $input_header_line;
    }
    print OUT $printed_input_header_line ;
    $output_header_line = $input_header_line;
  }
  
  if ($output_header_line ne $input_header_line) {
    print "Header line in file " . $input_file_name . " does not match first file:\n" . $input_file_names[0] . "\nSKIPPING FILE\n";
    continue;
  }
  
  while (my $line = <IN>) {
    if ($file_name_column) {
         $line = $base_file_name . ( ($format eq 'CSV') ? "," : "\t") . $line;
    }
    print OUT $line;
  }
  
  close IN;
}

#!/usr/bin/env perl -w

###
# Pod Documentation
###

=head1 NAME

batch_run.pl

=head1 SYNOPSIS

Usage: gene_gc_length.pl "[command]"

Run a command within each directory contained within the current directory. 
Enclose the command in quotes. All appearances of '#d' in the command will
be replaced with the name of the current directory.

=head1 DESCRIPTION

=over

=item B<-p> <pattern> 

Only execute in directories that match this regular expression.
Directories beginning with a period '.' or underscore '_' are always ignored.

=item B<-t> 

Test mode. Print commands that would have been run.

=item B<-c> 

Alternate way of passing the command to be run.

=back

=head1 AUTHOR

Jeffrey Barrick

=head1 COPYRIGHT

Copyright 2006.  All rights reserved.

=cut

###
# End Pod Documentation
###

use strict;

use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;

use Bio::SeqIO;

#Get options
use Getopt::Long;
use Pod::Usage;
my ($help, $man);
my ($input_file, $output_file);
pod2usage(1) if (scalar @ARGV == 0);
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input|i=s' => \$input_file,
	'output|o=s' => \$output_file,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

pod2usage(1) if (!$input_file);
pod2usage(1) if (!$output_file);

open OUT, ">$output_file";
print OUT join(",", ("gene", "length", "gccontent") ) . "\n";

my $in = Bio::SeqIO->new( -file   => "$input_file", -format => "GENBANK");
while (my $seq_object = $in->next_seq) {

  for my $feat_object ($seq_object->get_SeqFeatures) {
  
    my @line_list = ();
    if ($feat_object->primary_tag eq "gene") {

      push @line_list, $feat_object->get_tag_values('locus_tag');
      
      
      #Calculate gene length
      my $seq = $feat_object->seq()->seq;
      #print $seq . "\n";
      push @line_list, length $seq;
      
      
      my $gc = 0;
      #Calculate GC content
      my @seq_chars = split //, $seq;
      for my $c (@seq_chars) {
        if ( ("\U$c" eq "G") || ("\U$c" eq "C") ) {
          $gc++;
        }
      }
      
      $gc /= length $seq;
      push @line_list, $gc;
      
      print OUT join(",", @line_list) . "\n";
    }
  }
}

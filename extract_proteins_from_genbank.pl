#!/usr/bin/env perl -w

###
# Pod Documentation
###

=head1 NAME

extract_proteins_from_genbank.pl

=head1 SYNOPSIS

Usage: extract_proteins_from_genbank.pl -i input.gbk -o output.fna

Reads in a GenBank file and outputs a FASTA file of all protein sequences.

=head1 DESCRIPTION

=over

=item B<-i> <filename> 

Input GenBank file.

=item B<-o> <filename> 

Output FASTA file of protein sequences.

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
$output_file = "output.fna" if (!$output_file);

open OUT, ">$output_file";

my $in = Bio::SeqIO->new( -file   => "$input_file", -format => "GENBANK");
while (my $seq_object = $in->next_seq) {

  for my $feat_object ($seq_object->get_SeqFeatures) {
  
    if ($feat_object->has_tag('translation')) {
      
      my @locus_tags =  $feat_object->get_tag_values('locus_tag');
      my @translations =  $feat_object->get_tag_values('translation');
      print OUT ">$locus_tags[0]\n$translations[0]\n";
    }
  }
}

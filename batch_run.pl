#!/usr/bin/env perl

###
# Pod Documentation
###

=head1 NAME

batch_run.pl

=head1 SYNOPSIS

Usage: batch_run.pl "[command]"

Run a command within each directory contained within the current directory
(default more) or run a command for each file in a directory (if file mode used). 
Enclose the command in quotes. All appearances of '#d' in the command will be 
replaced with the name of the current directory/file. All appearances of '#e' 
will be replaced by the name with all extensions removed.

Type "perldoc batch_run.pl" or "batch_run.pl --man" for more detailed help!

=head1 DESCRIPTION

=over

=item B<-p,--pattern> <pattern> 

Only execute for directories/files that match this regular expression.
Names beginning with a period '.' or underscore '_' are always ignored.

=item B<-0,--file> 

File mode. Execute for each file in the current directory.

=item B<-c,--command> 

Alternate way of passing the command to be run as an option.

=item B<-t> 

Test mode. Print commands that would have been run. Don't execute them.

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
my ($command, $test, $pattern, $per_file);
#pod2usage(1) if (scalar @ARGV == 0);
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'command|c=s' => \$command,
	'test|t' => \$test,
	'pattern|p=s' => \$pattern,
	'file|0' => \$per_file,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

$command = "@ARGV";
print "Command $command\n";
print "TEST MODE -- no command executed!\n" if ($test);

my $current = `pwd`;
chomp $current;
print "Current directory: $current\n";

#Check out the directory that we are in for directories.
opendir(DIR, $current);
my @d = readdir DIR;
@d = grep !/^\./, @d; #Skip period prefixed
@d = grep !/^_/, @d; #Skip underscore prefixed
@d = grep /\Q$pattern\E/, @d if $pattern;

@d = grep {-d $_} @d if (!$per_file);

my $n = 0;
foreach my $d (@d)
{
	$d =~ s/\|/\\\|/g;
	
	$n++;

    my $temp_command = $command;
    #allow directory name as parameter in command: use '$d'
    $temp_command =~ s/\#d/$d/ge;
    #allow directory index as parameter in command: use '$n'
    $temp_command =~ s/\#n/$n/ge;
    
	if (!$per_file)
	{
   		print "cd $d; $temp_command\n";
		my $res = 0;
		$res = system "cd $d; $temp_command" if (!$test);
		#die "Command returned non-zero result code ($res). Exiting...\n" if ($res);
	}
	else
	{
		#allow file name without ending as parameter in command: use '$n'
		my $e = $d;
		$e =~ s/\..+//;
	    $temp_command =~ s/\#e/$e/ge;
		
		print "$temp_command\n";
		my $res = 0;
		$res = system "$temp_command" if (!$test);
	}
}


#!/usr/bin/perl 
####################################################
# parse PSIBLAST result file;
# INPUT: The psiblast output result.and the outputlist file. 
# OUPUT: A list of id and evalue, and score.
###################################################
use strict;

if (@ARGV != 2)
{
        die "INPUT: The psiblast output result.and the outputlist file..\n";
}

my $id_format = 30; # like AT1G18790.1  the formate of the ids of the sequence; suggested to be longer than the real id, just not too long to reach the score and evalue part
my $input_file = $ARGV[0];
my $output_file_2 = $ARGV[1];


open(READ, "<$input_file");   #this time find the hits;
open(WRITE, ">$output_file_2");
close(WRITE);  #to avoid the cannot remove error;
my $flag = "false";
while(<READ>){
	my $line = $_;
	if(substr($line, 0, 43) eq "Sequences producing significant alignments:"){
		$flag = "true";
		`rm $output_file_2`;
		
	}
	if (substr($line, 0, 1) eq ">"){
		$flag = "false"
	}
	if ($flag eq "true"){
		if(substr($line, 0, 43) ne "Sequences producing significant alignments:" && substr($line, 0, 40) ne "Sequences used in model and found again:" && length($line) > 2 && substr($line, 0, 65) ne "Sequences not found previously or not previously below threshold:"){
                                open(WRITE, ">>$output_file_2");
				if (substr($line, 0, 30)){
                                        print WRITE "$line";
                                }
				close(WRITE);
		}  
	}
}

#!/bin/perl -w
###########################################################
#used to run psi-blast on swiss prot and parse the results.
#Input:
#Ouput:
#
#Author: ZHENG WANG, June 10th, 2008.
###########################################################

use strict;


use lib "/rose/space1/CAFA2016/CAFA_website_2016/script/MIME-Lite-2.117/lib/";
use email;

use Cwd;
use MIME::Lite;
if(@ARGV != 8){
	die "Need 6 parameters: the path and file name of the sequence file containing all query sequences, the path of the psiblast executable program, the path of the directory where the single sequences are stored, the output folder path where the psiblast results are stored, the path of the database fasta file, the .dat database file for swissprot, and the LOG file path, and the original fasta file for the database;";
}
my $seq_file_all = $ARGV[0];
my $psi_blast_exe = $ARGV[1];
my $query_seq_dir = $ARGV[2]; #where to save the sequence files for query;
my $output_dir = $ARGV[3]; #where to save the psiblast results;
my $db_blast_dir = $ARGV[4];
my $db_content_dir = $ARGV[5];
my $log_file = $ARGV[6];
my $db_fasta_file = $ARGV[7]; 


my $start_time = time();


#######(1) Read the overal sequence file and seperate it into many single sequence files########### 
open(SEQALL, "<$seq_file_all");
open(LOG, ">$log_file");
##############Read db into array#######
print LOG "Reading the database into memory...\n";
open(DB, "<$db_content_dir");
my @db = <DB>;
close(DB);
print LOG "Reading the database complete.\n";
######################################
my $id;
my $seq;
my $tmpid;
my $tmpseq;
my $first_seq = "true";
my(@lineID);
while(<SEQALL>){

	my $line = $_;
	$line =~ s/\n//;
	if(length($line) <= 0){
		next;
	}
	$tmpid  = $id; # store the last seqid.
	$tmpseq = $seq; # store the last seq. 
	if($line =~ />/)
	{
		if($first_seq eq "false")
		{
			my $single_seq_path = $query_seq_dir."/".$id.".fasta";
			open(WRITE, ">$single_seq_path");
			print WRITE ">$id\n$seq\n";
			close(WRITE);
			print LOG "$id\n";
			$seq = "";
			
			my $blast_out = "$output_dir"."/"."$id".".out";
			if(-e $blast_out){
				print $id."existed\n";
				goto EXIST;
			}
			############run psiblast on the query single sequence on the database##########
			my $psi_blast_command = "$psi_blast_exe -i $single_seq_path -d $db_blast_dir -j 3 >$output_dir"."/"."$id".".out";	
			print $psi_blast_command."\n";
			print LOG "PSI-BLAST ".$psi_blast_command."\n";
			system("$psi_blast_command");

			############parse the psiblast results#########################################
		        my $list_file = "$output_dir"."/"."$id".".list";
			my $parse_psiblast_command = "perl /rose/space1/CAFA2016/CAFA_website_2016/swissprot_psiblast/parse_blast_result.pl $output_dir"."/"."$id".".out"." $list_file";
			print $parse_psiblast_command."\n";
			print LOG "PARSE ".$parse_psiblast_command."\n";
			system("$parse_psiblast_command");

			###########Get the top 1 blast hit and search swissprot .dat database #########
			###########to get the cross reference##########################################
			open(READLIST, "$list_file");
			my $ac_number;
			my $e_value;
			while(<READLIST>){
				my $line = $_;
				$ac_number = substr($line, 3, 6);
				my @items_2 = split(/\s+/, $line);
				my $len = scalar(@items_2);
				$e_value = $items_2[$len - 1];
				last;
			}
			close(READLIST);
			my $swiss_prot_part_output_file = $output_dir."/".$id.".swiss_prot_DR";
			if($ac_number){
				#print "$ac_number $e_value\n";
				#my $parse_swiss_prot_command = "perl /home/zwang/data/soybean_function/scripts/get_swiss_prot_parts.pl $db_content_dir $ac_number $swiss_prot_part_output_file DR"; #the database cross reference in swissprot begins with "DR";
				#print "$parse_swiss_prot_command\n";
				#print LOG "GET DR ".$parse_swiss_prot_command."\n";
				#`$parse_swiss_prot_command`;
				#############Parse the swiss-prot dataset##################
				my $ac_no = $ac_number;
				my $output_file = $swiss_prot_part_output_file;
				my $tag = "DR";
				open(WRITE2, ">$output_file");
				print WRITE2 "AC   $ac_no\n";
				my $found = "false";
				foreach my $line (@db){
        				if($found eq "true" && substr($line, 0, 2) eq "ID"){
                				$found = "false";
						last;
        				}
        				if($found eq "true" && substr($line, 0, 2) eq $tag){
                				print WRITE2 $line;
        				}
        				if(substr($line, 0, 2) eq "AC" && $line =~ /$ac_no/){
                				$found = "true";
        				}
				}
				close(WRITE2);
				##########END##############################################
				##############get the sequence with the $ac_number########
				#print "Finding sequence.\n";
				#print LOG "FINDING SEQUENCE\n";
				#my $seq;
				#my $found = "false";
				#open(FASTA, "<$db_fasta_file");
                                #while(<FASTA>){
				#	my $line = $_;
				#	if($found eq "true" && substr($line, 0, 1) eq ">"){
				#		$found = "false";
				#		last;
				#	}
				#	if($found eq "true"){
				#		$seq = $seq.$line;
				#		$seq =~ s/\n//g;
				#	}
				#	if(substr($line, 0, 7) eq ">$ac_number"){
				#		$found = "true";
				#	}
				#}
				#close(FASTA);
				open(WRITE, ">>$swiss_prot_part_output_file");
				#print WRITE "AA   $seq\n";
				print WRITE "EV   $e_value\n";
				close(WRITE);
			}
			else{
				print LOG "NO BLAST HIT\n";
				open(BLANK, ">$swiss_prot_part_output_file");
				close(BLANK);
			}

EXIST:			print LOG "\n";		
		}
		if($first_seq eq "true"){
			$first_seq = "false";
		}
		# $id = substr($line, 1, 6);
		@lineID = split(/\s+/,$line);
		$id = substr($lineID[0],1);
		print LOG $id;

	}
	else{
		$seq = $seq.$line;
		$seq =~ s/\n//;
	}
}
#################################################################################### 

### blast the last protein 
my $single_seq_path = $query_seq_dir."/".$id.".fasta";
open(WRITE, ">$single_seq_path");
print WRITE ">$id\n$seq\n";
close(WRITE);
print LOG "$id\n";
$seq = "";

my $blast_out = "$output_dir"."/"."$id".".out";
if(-e $blast_out){
	print $id."existed\n";
	goto EXIST;
}
############run psiblast on the query single sequence on the database##########
my $psi_blast_command = "$psi_blast_exe -i $single_seq_path -d $db_blast_dir -j 3 >$output_dir"."/"."$id".".out";	
print $psi_blast_command."\n";
print LOG "PSI-BLAST ".$psi_blast_command."\n";
system("$psi_blast_command");

############parse the psiblast results#########################################
my $list_file = "$output_dir"."/"."$id".".list";
my $parse_psiblast_command = "perl /rose/space1/CAFA2016/CAFA_website_2016/swissprot_psiblast/parse_blast_result.pl $output_dir"."/"."$id".".out"." $list_file";
print $parse_psiblast_command."\n";
print LOG "PARSE ".$parse_psiblast_command."\n";
system("$parse_psiblast_command");

###########Get the top 1 blast hit and search swissprot .dat database #########
###########to get the cross reference##########################################
open(READLIST, "$list_file");
my $ac_number;
my $e_value;
while(<READLIST>){
	my $line = $_;
	$ac_number = substr($line, 3, 6);
	my @items_2 = split(/\s+/, $line);
	my $len = scalar(@items_2);
	$e_value = $items_2[$len - 1];
	last;
}
close(READLIST);
my $swiss_prot_part_output_file = $output_dir."/".$id.".swiss_prot_DR";
if($ac_number){
	#############Parse the swiss-prot dataset##################
	my $ac_no = $ac_number;
	my $output_file = $swiss_prot_part_output_file;
	my $tag = "DR";
	open(WRITE2, ">$output_file");
	print WRITE2 "AC   $ac_no\n";
	my $found = "false";
	foreach my $line (@db){
			if($found eq "true" && substr($line, 0, 2) eq "ID"){
					$found = "false";
			last;
			}
			if($found eq "true" && substr($line, 0, 2) eq $tag){
					print WRITE2 $line;
			}
			if(substr($line, 0, 2) eq "AC" && $line =~ /$ac_no/){
					$found = "true";
			}
	}
	close(WRITE2);
	open(WRITE, ">>$swiss_prot_part_output_file");
	#print WRITE "AA   $seq\n";
	print WRITE "EV   $e_value\n";
	close(WRITE);
}
else{
	print LOG "NO BLAST HIT\n";
	open(BLANK, ">$swiss_prot_part_output_file");
	close(BLANK);
}
print LOG "\n";		


close(SEQALL);
close(LOG);
print "Done!\n";
#######OVER (1)#############

my $finish_time = time();
my $Blast_diff_hrs = ($finish_time -$start_time)/3600;
print "Blast finish! Total time: <$Blast_diff_hrs> hrs!\n"; 

my $from = "CAFA-Profunc <jh7x3\@mail.missouri.edu>";
my $to = "jh7x3\@mail.missouri.edu";
my $smtp = "localhost";
my $subject = "!! Blast finished";
my $body = "Attention: Blast finish! \nInput file: $seq_file_all \nThe output directory is : $output_dir. \nTotal running time: <$Blast_diff_hrs> hrs!";
####################
	my $subject1 = $subject." 1";
    my $msg1 = MIME::Lite->new(
    		From    => "$from",
    		To      => "$to",
    		Subject => "$subject1",
    		Cc     => 'miq44@mail.missouri.edu',
    		Type    => 'multipart/mixed'
    ) or die ("cannot create MIME object!");

    $msg1->attach(
        Type     =>'TEXT',
        Data     =>"$body"
	); 

    $msg1->send or print "Error sending email, MIME::Lite->send failed: $!\n";
####################



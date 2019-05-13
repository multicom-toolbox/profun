#! /usr/bin/perl -w
=pod
You may freely copy and distribute this document so long as the copyright is left intact. You may freely copy and post unaltered versions of this document in HTML and Postscript formats on a web site or ftp site. Lastly, if you do something injurious or stupid
because of this document, I don't want to know about it. Unless it's amusing.
=cut
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 use POSIX qw/ceil/; #this is for function ceil()



  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
	  print "The CAFA2 gives species for all the query sequences, this script just filter the blast result *.list, get the protein AC , and search uniprot database, get the taxid of each protein, and then filter the query protein in *.list file, ignore the protein which is not in the species!\n";

      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/24/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME addr_uniprot_dat dir_blast_result(*.list, the target name is T+species_ID+number, such as T96060000001) dir_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat \n\n");

	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";
 
 my($addr_uniprot)=$ARGV[0];
 my($dir_list)=$ARGV[1];      # all protein list file
 my($dir_output)=$ARGV[2];

 -s $addr_uniprot || die "Cannot open $addr_uniprot!\n";
 -s $dir_list || die "Cannot open $dir_list!\n";
 -s $dir_output || system("mkdir $dir_output");

 
 my($key,$value,$IN,$OUT,$path_input,$target_name,$path_output,$line,$file,$i,$j,$evalue,$network_type);
 my(@tem_split,@tem_2,@files);


 $| = 1;
 print "Loading uniprot database and get the species for each protein ... ";
 #sleep(1);
 #############################################################################################
 #
 #                1. Load uniprot database
 #
 #############################################################################################
 my(%ID2spe)=();              # the protein ID mapped to species

 
 my($index_ID,$index_AC,$key_spe,$i1,$i2);
 my(@v_AC,@species_ID,@tem_for_spe);
 my(%tem_spe);

 @species_ID=();              # store all species ID for this 
 @v_AC=();
 $index_ID=0;
 $index_AC=0;

 $IN=new FileHandle "$addr_uniprot";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
     @tem_split=split(/\s+/,$line);
	 if(substr($line,0,2) eq "//")
	 {# new line comes
         if(@species_ID == 0)
		 {
			 print "Warning, check this @v_AC, cannot find species for this protein record, check!\n";
			 next;
		 }

         for($i=0;$i<@v_AC;$i++)
		 {# build relationship of AC and species
			 if(not exists $ID2spe{$v_AC[$i]})
			 {
				 $key_spe=$species_ID[0];
				 for($i1=1;$i1<@species_ID;$i1++)
				 {
					 $key_spe.="|".$species_ID[$i1];
				 }
				 $ID2spe{$v_AC[$i]}=$key_spe;
			 }
			 else
			 {
				 print "!!!!!! this protein AC $v_AC[$i] mapped to more than one species ID !!!!\n";
                 @tem_for_spe=split(/\|/,$ID2spe{$v_AC[$i]});
				 %tem_spe=();
				 for($i1=0;$i1<@tem_for_spe;$i1++)
				 {
					 if(not exists $tem_spe{$tem_for_spe[$i1]})
					 {
						$tem_spe{$tem_for_spe[$i1]}=1;
					 }
				 }
				 for($i1=0;$i1<@species_ID;$i1++)
				 {
					 if(not exists $tem_spe{$species_ID[$i1]})
					 {
						$tem_spe{$species_ID[$i1]}=1;
					 }
				 }
				 $j=0;
				 foreach $value (keys %tem_spe)
				 {
					 if($j==0)
					 {
						 $ID2spe{$v_AC[$i]}=$value;
						 $j=1;
					 }
					 else
					 {
						 $ID2spe{$v_AC[$i]}.="|".$value;
					 }
				 }
			 }# end else


		 }


		 @species_ID=();              # store all species ID for this 
		 @v_AC=();
		 $index_ID=0;
		 $index_AC=0;
         
	     next;
	 }
     if($tem_split[0] eq "AC")
	 {
		 for($i=1;$i<@tem_split;$i++)
		 {
			 @tem_2=split(/\;/,$tem_split[$i]);
			 $v_AC[$index_AC]=$tem_2[0];
			 $index_AC++;
		 }
	 }
	 elsif(($tem_split[0] eq "OX"))           # here we only consider OX
	 {# get the species, such as line is : OX   NCBI_TaxID=654924;
		 
		 @tem_for_spe=split(/\;/,$tem_split[1]);
		 @tem_for_spe=split(/\=/,$tem_for_spe[0]);

		 $species_ID[$index_ID] = $tem_for_spe[1];
		 $index_ID++;
	 }


 }

 $IN->close();
 print "Done!\n";
 print "Processing protein lists ...";
 #############################################################################################
 #
 #                2. Load protein lists
 #
 #############################################################################################
 
 my(@all_ids);
 $all_ids[0] = "186497";
 $all_ids[1] = "243232";
 $all_ids[2] = "273057";
 $all_ids[3] = "309800";
 $all_ids[4] = "436308";
 $all_ids[5] = "453591";
 $all_ids[6] = "478009";
 $all_ids[7] = "160488";
 $all_ids[8] = "170187";
 $all_ids[9] = "208964";

 $all_ids[10] = "223283";
 $all_ids[11] = "224308";
 $all_ids[12] = "243273";
 $all_ids[13] = "321314";
 $all_ids[14] = "83333";
 $all_ids[15] = "85962";
 $all_ids[16] = "99287";
 $all_ids[17] = "10090";
 $all_ids[18] = "10116";
 $all_ids[19] = "284812";
 $all_ids[20] = "3702";
 $all_ids[21] = "44689";
 $all_ids[22] = "559292";
 $all_ids[23] = "7227";
 $all_ids[24] = "7955";
 $all_ids[25] = "8355";
 $all_ids[26] = "9606";
 my(%all_possible_spe)=();     # get all possible species
 for($i=0;$i<@all_ids;$i++)
 {
	 $all_possible_spe{$all_ids[$i]}=1;
 }
 my($spe_id);


 opendir(DIR,"$dir_list");
 @files=readdir(DIR);
 foreach $file (@files)
 {
	 if($file eq "." || $file eq "..")
	 {
		 next;
	 }
	 @tem_split=split(/\./,$file);
	 if($tem_split[@tem_split-1] ne "list")
	 {
		 print "skip file $file, the name is not *.list\n";
		 next;
	 }
     ####### get the taxid of this file #######
	 $i=length($tem_split[0])-8;
	 $spe_id=substr($tem_split[0],1,$i);         # get the species name
	 if(not exists $all_possible_spe{$spe_id})
	 {
		 print "Fails, check the file $file, the species ID is $spe_id, not possible!\n";
		 exit(0);
	 }
     ##########################################

	 $path_input=$dir_list."/".$file;
	 $path_output=$dir_output."/".$file;      # this is the filtered output file
	 $OUT = new FileHandle ">$path_output";
	 $IN=new FileHandle "$path_input";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;
         @tem_split=split(/\|/,$line);
		 if(not exists $ID2spe{$tem_split[1]})
		 {
			 print "Warning, the $tem_split[1] is not in any species!\n";
			 next;
		 }
		 else
		 {
			 if($ID2spe{$tem_split[1]} eq $spe_id)
			 {
				 print $OUT $line."\n";
			 }
			 else
			 {
				 print "Skip $line, $spe_id is not the same as $ID2spe{$tem_split[1]}!\n";
				 next;
			 }
		 }
	 }
	 $IN->close();
	 $OUT->close();
 }
 
 

  my $endtime = localtime();
  print "\nThe time ended at : $endtime.\n";

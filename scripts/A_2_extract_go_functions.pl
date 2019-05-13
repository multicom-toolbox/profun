##########################################################################################################
#    Function about getting all GO terms of contact gene pairs,                     check GO in clster.  #
#                                  output: addr_output                                                   #
#										Renzhi Cao  													 #            
#																										 #            
#									    11/25/2012														 #            
#																										 #            
#																										 #            
#									Revised at 4/12/2013												 #            
#																										 #            
##########################################################################################################
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
 sub add_sequence($);             # add a sequence to a sequence vector.

  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args

	 
      print("This program will select protein pairs of irefindex database processed network(two protein IDs are included in the input network), also check uniprot to make sure all protein has GO functions. output the GO functions of these protein pairs. Output format is GO:***** GO:****** for each line, and #ENDOFGENE,#ENDOFCLUSTER.\n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 9/26/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME addr_swiss_prot_dat dir_input_ppis(*.network file, each file has format like:proteinID1 proteinID2 taxid1 taxid2) dir_output!\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../data/converted_network_with_taxid ../result/0_all_go_functions_for_ppi \n\n");

	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";
 my($addr_swiss)=$ARGV[0];
 my($dir_input)=$ARGV[1];

 my($dir_output)=$ARGV[2];

 -s $dir_input || die "No input file $dir_input!\n";
 -s $dir_output || system("mkdir $dir_output");

 my($key,$file,$value,$IN,$OUT,$line,$i,$j);
 my(@files,@tem_split,@tem_2);
 my(%hash_all_proteinID)=();                   # store all protein IDs used in the ppi

 my(%hash_sequence_of_protein)=();           # store all sequence of the protein
 
    my(%protein_GO)=();
    my(@tem_2);
	my ($to_search)=1; # 1 means it need to search
	my ($hit_protein) = 0; # to see whether this protein is what we want
	my($hit_go)=0;

    
	my($hit_seq)=0;
	my($seq_vector)="NULL";

	my($GO_vector)=":";
	my($saved_protein);
    my($iii);
    my(@array_contact_pairs)=();            # store gene pairs
    my(@go_1,@go_2,@filtered_go_1,@filtered_go_2);
    my($index_array)=0;
    my($input,$addr_output);

    my($random_range)=$index_array;
    my($random_total)=0;              # count the random number of gene pairs we select
    my($random_selected,$random_selected_pair,$genei,$genej,$k1,$k2);

######### 0 . Get each ppi network, and start doing the same thing ##############
 opendir(DIR,"$dir_input");
 @files=readdir(DIR);
foreach $file (@files) 
{
	if($file eq "." || $file eq "..")
	{
		next;
	}
        @tem_split=split(/\./,$file);
        if($tem_split[@tem_split-1] ne "network")
        {
          print "Warning !!! We skip $file, the name is not *.network!\n";
          next;
        }
	$input = $dir_input."/".$file;
	$addr_output = $dir_output."/".$file;
    %hash_all_proteinID=();
	%hash_sequence_of_protein=();





    print "Processing $input ...\n";
 #################################################################
 #
 #           1. Get all proteinID
 #
 #################################################################
 $IN=new FileHandle "$input";
 defined($IN) || die "Cannot open input $input!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<2)
	 {# this is a ppi network file, so has protein1 protein2
		 next;
	 }
	 if(!exists $hash_all_proteinID{$tem_split[0]})
	 {
		 $hash_all_proteinID{$tem_split[0]}=1;
	 }
     if(!exists $hash_all_proteinID{$tem_split[1]})
	 {
		 $hash_all_proteinID{$tem_split[1]}=1;
	 }

	 if(!exists $hash_sequence_of_protein{$tem_split[0]})
	 {
		 $hash_sequence_of_protein{$tem_split[0]}="NULL";
	 }
     if(!exists $hash_sequence_of_protein{$tem_split[1]})
	 {
		 $hash_sequence_of_protein{$tem_split[1]}="NULL";
	 }


 }
 $IN->close();
 print "All protein ID has been loaded!\n";
 #################################################################
 #
 #           2. Search swiss prot to get all GO functions for all protein AC numbers
 #
 ################################################################# 
    %protein_GO=();

	$to_search=1; # 1 means it need to search
	$hit_protein = 0; # to see whether this protein is what we want
	$hit_go=0;

    
	$hit_seq=0;
	$seq_vector="NULL";

	$GO_vector=":";
	
    $IN = new FileHandle "$addr_swiss";
    if (! defined($IN)) 
    {
        print("Can't open spec file $addr_swiss: $!\n");
       return;
    }
	while ( defined($line = <$IN>))
	{
	   $line =~ s/\s+$//;
       chomp($line); # strip off trailing '\n' (newline)
       if($line eq "//")
	   {
############### save the GO function for each AC#####################
           if($hit_protein==1 && $hit_go==1)
		   {
			   if(! exists $protein_GO{$saved_protein})
			   {
				   $protein_GO{$saved_protein}=$GO_vector;
			   }
			   ########## save sequence for each AC #########
			   if(! exists $hash_sequence_of_protein{$saved_protein})
			   {
				   print "The protein is not in the hash list???? Error!\n";
			   }
			   $hash_sequence_of_protein{$saved_protein}=$seq_vector;
		   }
		   
############### save the GO function for each AC#####################
           $hit_seq=0;
		   $seq_vector="NULL"; 

		   $to_search = 1;
		   $hit_protein=0;
		   $hit_go=0;
		   $GO_vector=":";
		   next;
	   }
	   if($to_search == 1)
	   {#we will search
		   @tem_split = split(/\s+/,$line);
		   if($tem_split[1] eq "GO;")
		   {
			   @tem_split=split(/\;/,$line);
			   if($GO_vector eq ":")
			   {
				   $i=1;
				   @tem_2=split(/\s+/,$tem_split[$i]);
				   $GO_vector=$tem_2[1];
				   for($i=4;$i<@tem_split;$i+=3)
				   {
					   @tem_2=split(/\s+/,$tem_split[$i]);
                       $GO_vector=$GO_vector."|".$tem_2[1];
				   }
			   }
			   else
			   {
				   for($i=1;$i<@tem_split;$i+=3)
				   {
					   @tem_2=split(/\s+/,$tem_split[$i]);
                       $GO_vector=$GO_vector."|".$tem_2[1];
				   }
			   }
               $hit_go=1;
		   }
		   elsif($tem_split[0] eq "AC")
		   {  # get the AC number line, now try to search on this line
			   
			   for($iii = 1; $iii < @tem_split; $iii++)
			   {
				  @tem_2=split(/\;/,$tem_split[$iii]);
			      $saved_protein=$tem_2[0];
			      if(exists $hash_all_proteinID{$saved_protein})
			      {# we already hit one protein ID . 
                      $hit_protein=1;
					  last;
			      }
			   }
               
			   if($iii == @tem_split)
			   {# we don't find the protein AC number we needed
				   $to_search=0;
			   }
		   }
		   elsif($tem_split[0] eq "SQ")
		   {# this is for getting the sequence
			   $hit_seq=1;
		   }
		   else
		   {
			   if($hit_seq == 1)
			   {# this is for sequence 
				   if($seq_vector eq "NULL")
				   {
					   $seq_vector=add_sequence($line);
				   }
				   else
				   {
					   $seq_vector.=add_sequence($line);      # add the sequence to the sequence vector, if seq_vector is "NULL", then initialize it, otherwise, add the new line.

				   }
				   
			   }
		   }
	   }
	}#end while
    $IN->close();
print "Swiss prot database checked!\n";
 #################################################################
 #
 #           3. Read protein pairs file again, check the protein pairs which have function,output the GO term
 #
 #################################################################
  if(-e $addr_output)
  {
    print "the result file | : $addr_output  ...Exists!\n"; 
  }
  else 
  { 
    open (File, "&gt;$addr_output");
    chmod (0777, $addr_output); 
    close (File);
  }      
  $OUT = new FileHandle "> $addr_output";
  defined($OUT) || die "Cannot open output file $addr_output!\n";

 @array_contact_pairs=();            # store protein pairs

 $index_array=0;
 $IN=new FileHandle "$input";
 defined($IN) || die "Cannot open input $input!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<2)
	 {# this is a ppi network file, so has protein1 protein2
		 next;
	 }
     $array_contact_pairs[$index_array]=$tem_split[0]."|".$tem_split[1];
	 $index_array++;
 }
 $IN->close();



 $random_range=$index_array;


 for($i=0;$i<$index_array;$i++)
 {#
	 
     $random_selected_pair=$array_contact_pairs[$i];    # the protein pair selected
	 @tem_split=split(/\|/,$random_selected_pair);
	 $genei=$tem_split[0];
	 $genej=$tem_split[1];
     if((exists $protein_GO{$genei}) && (exists $protein_GO{$genej}))
	 {# this gene pair has GO function
		 # output the GO function for this two genes
         @go_1=split(/\|/,$protein_GO{$genei});
         @go_2=split(/\|/,$protein_GO{$genej});
		 
         @filtered_go_1=();
		 @filtered_go_2=();
		 $k2=0;
		 for($k1=0;$k1<@go_1;$k1++)
		 {
			 #if(substr($go_1[$k1],0,6) eq "GO:000")
			 #{
				 $filtered_go_1[$k2]=$go_1[$k1];
				 $k2++;
			 #}
		 }
		 if($k2==0)
		 {# no function we need
			 next;
		 }
		 $k2=0;
		 for($k1=0;$k1<@go_2;$k1++)
		 {
			 #if(substr($go_2[$k1],0,6) eq "GO:000")
			 #{
				 $filtered_go_2[$k2]=$go_2[$k1];
				 $k2++;
			 #}
		 }
		 if($k2==0)
		 {# no function we need
			 next;
		 }
		 ###########################################################################################
         ##### output the pairwise GO function ######
		 for($k1=0;$k1<@filtered_go_1;$k1++)
		 {
			 for($k2=0;$k2<@filtered_go_2;$k2++)
			 {
				 print $OUT $filtered_go_1[$k1]."\t".$filtered_go_2[$k2]."\t"."1"."\n";
			 }
		 }
         print $OUT "#ENDOFGENE\t#ENDOFGENE $genei|$hash_sequence_of_protein{$genei}|$genej|$hash_sequence_of_protein{$genej}\n";
		 $random_total++;
	 }
	 else
	 {# at least one of them has no GO function, skip
		 next;
	 }
 }
 $OUT->close();

}

 sub add_sequence($)
 {# extract the sequence from a line. 
	 my($line)=@_;
	 my(@tem_split)=split(/\s+/,$line);
	 if($tem_split[0] ne "")
	 {
		 print "Check this line : $line\n";
		 print "Maybe the sequence is wrong??? \n";
		 exit(0);
	 }
	 my($seq);
	 $seq=$tem_split[1];
	 my($i);
	 for($i=2;$i<@tem_split;$i++)
	 {
		 $seq.=$tem_split[$i];
	 }
	 $seq=uc($seq);                # convert to uppercase
	 return $seq;
 }

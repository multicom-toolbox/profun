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
 sub get_one_step_neighbour($);     # get one step neighbor from the input line.

  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
      print("This program will check the mapped protein ID/gene ID, and check the network to get one step neighbors for that network. \n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/6/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir_one_step_neighbours(*.one_step_gene_neighbours or *.one_step_protein_neighbours) protein_list(proteinID evalue_target_name species mapped_protein/gene) addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME ../result/3_one_step_protein_neighbours ../result/F_2_protein_lists_with_network_type ../result/F_3_protein_lists_with_network_and_one_step_neighbour \n\n");

	  exit(1) ;
    }
  
   my($dir_neighbour)=$ARGV[0];
   my($addr_list)=$ARGV[1];
   my($output)=$ARGV[2];
   
   -s $dir_neighbour || die "Cannot open the directory of one step protein/gene neighbors Check script A_5_get_one_step_gene_neighbours_from_hi_c_gene_network.pl and A_5_get_one_step_protein_neighbours_from_ppi.pl!\n";
   -s $addr_list || die "Cannot open protein lists file, check F_1_get_all_proteins_from_blast_result.pl\n";

   my($IN,$OUT,$line,$file,$key,$n_function,$i,$path_input,$type);
   my(@files,@tem_split,@tem_2);
   my(%hash_neighbours)=();
##############################################################################
#
#
#               1. Load all protein neighbors or gene neighbors 
#
#
##############################################################################
   opendir(DIR,"$dir_neighbour");
   @files=readdir(DIR);
   foreach $file (@files)
   {
       if($file eq "." || $file eq "..")
       {
           next;
       }
       @tem_split=split(/\./,$file);
       if(@tem_split<2)
       {
           print "warning! We need the file name *.one_step_gene_neighbours or *.one_step_protein_neighbours, skip $file!\n";
           next;
       }
       $path_input=$dir_neighbour."/".$file;
       $type=$tem_split[0];          # the network type
       #### now read the one step neighbors and store that in the hash table ####
       $IN=new FileHandle "$path_input";
       while(defined($line=<$IN>))
       {
           chomp($line);
           @tem_split=split(/\|/,$line);
           if(@tem_split<1)
           {
              next;
           } 
           $key=$type."|".$tem_split[0];         # network type and the ID, such as 3_NORMAL_gene|GeneID:64168
           if(not exists $hash_neighbours{$key})
           {
               @tem_2=split(/\s+/,$tem_split[1]);
               if(@tem_2<1)
               {
                   print "Special warning, there is no neighbor for this ID $tem_split[0]\n";
                   next;
               }
               $n_function=$tem_2[0];
               for($i=1;$i<@tem_2;$i++) 
               {
                   $n_function.="_".$tem_2[$i];
               }
               $hash_neighbours{$key}=$n_function;
           }
           else
           {
               print "Very strange thing, this $key has multiple places , $type, check !\n";
           }
       }
       $IN->close();
      
   }
 
##############################################################################
#
#
#               2. process the input protein lists, add the one step protein/gene neighbors to it 
#
#
##############################################################################
   $OUT = new FileHandle ">$output"; 

   my($neighbours);
   $IN=new FileHandle "$addr_list";
   while(defined($line=<$IN>))
   {
      chomp($line);
      $line=~s/\s+$//;
      @tem_split=split(/\s+/,$line);
      if(@tem_split<4)
      {
          print "Not correct format(proteinID evalue_target_name species mapped_protein/gene). check $line\n";
          next;
      }
      $neighbours=get_one_step_neighbour($line);
      print $OUT $line."\t".$neighbours."\n";
   }
   $IN->close();

   $OUT->close();


sub get_one_step_neighbour($)
{
   my($input_line)=@_;
   my(@mapped_species)=();    # the mapped species for the protein
   my(@mapped_IDs)=();        # the mapped proteinID/geneID for the query protein
   my($key,$mapped_neighbour);
   @tem_split=split(/\s+/,$input_line);
   if($tem_split[2] eq "NULL" || $tem_split[3] eq "NULL")
   {
       $mapped_neighbour="NULL";
       return $mapped_neighbour;
   }

   @mapped_species=split(/\;/,$tem_split[2]);
   @mapped_IDs=split(/\;/,$tem_split[3]);              # the mapped_IDs[$i] could have multiple ID

   my(@spe_ID);

   my($len1,$i,$len2);
   $len1=scalar(@mapped_species);
   $len2=scalar(@mapped_IDs);
   if($len1 != $len2)
   {
      print "Wrong format, check $input_line! in the protein list, the number of species mapped should be identical to the number of protein/gene IDs\n"; 
      exit(0);
   }
   
   for($i=0;$i<@mapped_IDs;$i++)
   {
	   @spe_ID=split(/\|/,$mapped_IDs[$i]);
	   if(@spe_ID>1)
	   {
		   print "This protein mapped to multi proteins, check $mapped_IDs[$i], the input line is $input_line! we only use one of them, because they are actually the same record\n";
		   $mapped_IDs[$i]=$spe_ID[0];
	   }
   }

   $key=$mapped_species[0]."|".$mapped_IDs[0];       # get the species and IDs
   if(not exists $hash_neighbours{$key})
   {
       print "Warning,  there is no neighbor for $key, could be wrong, check the reason!\n";
       $mapped_neighbour="NULL";
   }
   else
   {
      $mapped_neighbour=$hash_neighbours{$key};
   }
   for($i=1;$i<@mapped_IDs;$i++)
   {
       $key=$mapped_species[$i]."|".$mapped_IDs[$i];
       if(not exists $hash_neighbours{$key})
       {
           print "Could be wrong, not exists neighbour for $key!\n";
           $mapped_neighbour.="|NULL";
       }
       else
       {
            $mapped_neighbour.="|".$hash_neighbours{$key};
       }
   }
   return $mapped_neighbour;
}


#For all mapped protein ID/ gene ID, we check the network, and get one step neighbours for that.
















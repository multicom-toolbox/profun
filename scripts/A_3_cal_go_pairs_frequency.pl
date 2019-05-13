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


   if (@ARGV != 2)
    { # @ARGV used in scalar context = number of args
          print "I update this script, to make it faster !!!!\n";
          print "For the hi-C network, check B_3_cal_go_pairs_frequency.pl\n";
	  print "This script generates the go pair frequency from an input go pair file, you can use the script in 4_10_script/0_random_select_gene_pairs_more_than_threshold.pl to get the go pairs from the gene network. !\n";
	  print "\n************** Renzhi Cao *******************\n";
	  print "Input:\n";
	  print "1. Input file with the go pairs \n";
	  print "2. Address of output. (All the go pairs frequency existing in the gene pairs with high contact)\n";
      print "\n************** Renzhi Cao *******************\n";
	  print "For example:\n";
          #print "For CAFA2:\n";
          print "perl $0 ../result/0_all_go_functions_for_ppi/other.network ../result/1_go_frequency/other.go_frequency\n";

print "perl B_3_cal_go_pairs_frequency.pl ../result/B_2_gene_network_go_function_pairs/all_go_pairs_for_3_NORMAL ../result/B_3_go_frequency/3_NORMAL.go_frequency\n";

      exit(0);
	}


############### the copy right for each script #############################
print  "REM       \n";
print  "REM          *********************************************************************************\n";
print  "REM          *                                                                               *\n";
print  "REM          *                                                                               *\n";
print  "REM          *      Developed By :   Renzhi Cao (rcrg4\@mail.missouri.edu)                   *\n";
print  "REM          *      Copyright    :   Dr. Jianlin Cheng's BDML Laboratory                     *\n";
print  "REM          *      Release Date :   June 21, 2013                                          *\n";
print  "REM          *      Vesion       :   1.0                                                     *\n";
print  "REM          *                                                                               *\n";
print  "REM          *********************************************************************************\n";
print  "REM       \n";
############################################################################
    my($input)=$ARGV[0];      # input of go pairs
	my($output)=$ARGV[1];     #  output file

  my(%hash_GO)=();            # GO partner frequency
#  my(%hash_all_go)=();        # key is GO, value is the frequency of other GO and the frequency of other GO, which is in contact with this GO.
  
  ######### read the go pair, first get the GO pairs frequency ###########
  my($IN,$OUT,$line,$key,$key2,$i);

  
  my($tmp_folder)=$output."_tmp";
  #### create a tmp folder to store the frequency of GO partner for each GO item
  if(!-s $tmp_folder)
  {
     system("mkdir $tmp_folder");
  }
  else
  {
      system("rm -R $tmp_folder");
      system("mkdir $tmp_folder");
  }
  
  my(@tem_split,@tem222);
  my($GO_path);
  $IN = new FileHandle "$input";
  while(defined($line=<$IN>))
  {
	  chomp($line);
	  $line=~s/\s+$//;
	  @tem_split=split(/\s+/,$line);
	  if(@tem_split<2)
	  {
		  next;
	  }
	  if($tem_split[0] eq "#ENDOFGENE")
	  {
		  next;
	  }
######################################################## 
#           Store the frequency of GO partner for the GO, use the GO term as the file name 
          $GO_path=$tmp_folder."/".$tem_split[0];
          $OUT = new FileHandle ">>$GO_path";
          print $OUT $tem_split[1]."\t".$tem_split[2]."\n";
          $OUT->close();
          $GO_path=$tmp_folder."/".$tem_split[1];
          $OUT = new FileHandle ">>$GO_path";
          print $OUT $tem_split[0]."\t".$tem_split[2]."\n";
          $OUT->close();

#########################################################
  }
  $IN->close();
######   now collect each GO file, calculate the frequency and output the result ####
  my(@files);
  my($file,$tag_1);
  $OUT = new FileHandle ">$output";
  my(@array)=();
  my(@fre)=();
  my($index)=0;
  my($count_number,$max);
  $count_number=0;
  opendir(DIR,"$tmp_folder");
  @files=readdir(DIR);

  my($total_number)=scalar(@files);
  foreach $file (@files)
  {
      if($file eq "." || $file eq "..")
      {
           next;
      }

      $count_number++;
      print "$count_number/$total_number...";
          

      $GO_path=$tmp_folder."/".$file;
      %hash_GO=();
      $IN = new FileHandle "$GO_path";
      while(defined($line= <$IN>))
      {
         chomp($line);
         @tem_split=split(/\s+/,$line);
         if(@tem_split<2)
         {
             next;
         }
         if(not exists $hash_GO{$tem_split[0]})
         {
             $hash_GO{$tem_split[0]}=$tem_split[1];
         }
         else
         {
             $hash_GO{$tem_split[0]}+=$tem_split[1];
         }
      }
      $IN->close();
      ###### rank all GO by their frequency and output the result ######
      @array=();
      @fre=(); 
      $index=0;
      foreach $key (keys %hash_GO)
      {
          for($i=$index;$i>0;$i--)
          {
               if($hash_GO{$key}>$fre[$i-1])
               {
                   $array[$i]=$array[$i-1];
                   $fre[$i]=$fre[$i-1];
               }
               else
               {
                   last;
               }
          }
          $array[$i]=$key;
          $fre[$i]=$hash_GO{$key};
          $index++;
      }

      print $OUT $file."\t";
      if($index!=0)
      {
          print $OUT $array[0]."_".$fre[0];
      }
      else
      {
          print "Strange, this $file has no partner! check this file!!!!\n";
          exit(0);
          #next;
      }
      for($i=1;$i<$index;$i++)
      {
          print $OUT "|".$array[$i]."_".$fre[$i];
      }
      print $OUT "\n";  
  }

  $OUT->close();

  system("rm -R $tmp_folder");

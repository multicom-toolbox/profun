#! /usr/bin/perl -w
=pod
You may freely copy and distribute this document so long as the copyright is left intact. You may freely copy and post unaltered versions of this document in
 HTML and Postscript formats on a web site or ftp site. Lastly, if you do something injurious or stupid
because of this document, I don't want to know about it. Unless it's amusing.
=cut
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 use POSIX qw/ceil/; #this is for function ceil()
 sub get_num_columns($);
 sub get_num_rows($);
 sub get_all_headers($);
 sub get_max_len($);       # get the maximum length of target name
  if (@ARGV < 2)
    { # @ARGV used in scalar context = number of args
          print("This program try to create a R script (TMP.r) in the output file to draw svg file for the input matrix\n");
		  print("***********************************************************\n\n");
          print("Renzhi Cao, 3/7/2012\n\n");
          print("***********************************************************\n");
          print("You should execute the perl program like this: perl $PROGRAM_NAME addr_input(the file should be four columns (T859620000003   GO:0046872      0.75    molecular_function) for each type, without header) addr_output !\n");
		  print "Please check inside the script, to change the input title or any other sub titles you need![xlab_titles ylab_titles X_axis_sub_labels1 sub_labels2 sub_labels3 ...]\n";
          print("\n\n***********************************************************\n");
          print("Examples:\n");
		  print("perl $PROGRAM_NAME ../test/R_all_output_08110907_v2/S1_MHS_net_seq_score ../test/f1.svg\n\n");
          exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";



  my($addr_input)=$ARGV[0];

  my($addr_output)=$ARGV[1];


  my($X_title)="NULL";
  my($Y_title)="NULL";
  if(@ARGV>3)
  {
          $X_title=$ARGV[2];
          $Y_title=$ARGV[3];
  }
  my(@sub_xlab)=();
  my($i_x,$index_i);
  if(@ARGV>5)
  {
          $i_x=0;
          for($index_i=4;$index_i<@ARGV;$index_i++)
          {
                  $sub_xlab[$i_x++]=$ARGV[$index_i];
          }
  }


  my($tmp_output)=$addr_output.".TMP.r";
#############################################################
  my(@tem_split,$line,$i,$j,$k);
  my($IN,$OUT);
  my($path,$output_path,$num_columns,$num_rows,$max_len);
  my($file);
  my(@files);
  my(@all_colors)=("blue","red","forestgreen","darkgreen","brown1","cyan","darkgoldenrod","darkseagreen","gold","lawgreen");        # we only use 10 colors here, revise here
  my(@all_headers)=();
######################## open the input dir #################

##########do something to the file###################
        $path=$addr_input;
        $output_path=$addr_output;
        $max_len = get_max_len($path);                       # get the max length of the target name
        $num_columns=get_num_columns($path);                 # get the number of columns from input file
        if($num_columns != 4)
		{
			print "Warning, there should be only four columns! check the input file\n";
			
		}

        print "We catch file : $path\n";
#################### output for r script ################################
        open (File, "&gt;$tmp_output");
        chmod (0777, $tmp_output);
        close (File);
        $OUT = new FileHandle "> $tmp_output";
        if (! defined($OUT) )
        {
           croak "Unable to open output file: $tmp_output. Bye-bye.";
           exit(1);
        }
############ write a r script ############################
    ########### set the margin ###########
        #print $OUT "par(mar=c(5,116,4,2)+0.1,mgp=c(15,1,0))\n";


        print $OUT "x=read.table(\"$path\", header=F)\n";
        print $OUT "names(x) <- c(\"target\",\"GO\",\"score\",\"category\")\n";
		print $OUT "x\$category=gsub(\"molecular_function\",\"MF\",x\$category)\n";
		print $OUT "x\$category=gsub(\"biological_process\",\"BP\",x\$category)\n";
		print $OUT "x\$category=gsub(\"cellular_component\",\"CC\",x\$category)\n";
        print $OUT "targetname <- x\$target\n";
		print $OUT "category <- x\$category\n";
		print $OUT "fun <- data.frame(targetname,category)\n";
        print $OUT "svg(\"$output_path\", height=10, width=10)\n";



    ########## change the mar and mpg #################
    #print $OUT "par(mar=c(5.1,7.1,4.1,2.1),mgp=c(4.5,1,0))\n";
        print $OUT "par(mar=c(15,4.1,4.1,2.1))\n";
        ###### set the xlab vertical ########
        print $OUT "par(las=2)\n";
    ############# plot the main dot lines ############
        print $OUT "spineplot(category~targetname,data = fun, breaks = 3, col = c(\"red\",\"green\",\"blue\"), xlab=\"\")\n";
        print $OUT "par(las=1)\n";
        print $OUT "mtext(\"target name\",side=1,line=$max_len/2+3)\n";
        print $OUT "dev.off()\n";
        $OUT->close();
        system("R --vanilla < $tmp_output");
  
  system("rm $tmp_output");
  ##########################################################################################
  my $endtime = localtime();
  print  "\n$PROGRAM_NAME ended at : $endtime.\n";

  ##########################################################################################################
#      Function about transform the two dimentional coordination to linear array coordination            #
#                                                                                                                                                            #
#                                                                               Renzhi Cao                                                                   #
#                                                                                                                                                            #
#                                                                           12/27/2011                                                                       #
#                                                                                                                                                            #
##########################################################################################################


 sub get_num_columns($)
 {
         my($input)=@_;
         my($IN,$line);
         my(@tem_split);
         $IN= new FileHandle "$input";
         if(defined($line=<$IN>))
         {
                 chomp($line);
                 @tem_split=split(/\|/,$line);
         }
         $IN->close();
         my($i)=scalar(@tem_split);
         return $i;
 }
 sub get_all_headers($)
 {
         my($input)=@_;
         my($IN,$line);
         my(@tem_split);
         $IN= new FileHandle "$input";
         if(defined($line=<$IN>))
         {
                 chomp($line);
                 @tem_split=split(/\|/,$line);
         }
         $IN->close();

         return @tem_split;
 }
  sub get_num_rows($)
 {
         my($input)=@_;
         my($IN,$line);
         my(@tem_split);
         my($i)=0;
         $IN= new FileHandle "$input";
         while(defined($line=<$IN>))
         {
                 chomp($line);
                 $i++;

         }
         $IN->close();
         $i--;
         return $i;
 }
 sub get_max_len($)
 {
      my($input)=@_;
      my($IN,$line);
      my(@tem);
      my($len)=0;
      $IN= new FileHandle "$input";
      while(defined($line=<$IN>))
      {
          chomp($line);
          @tem=split(/\s+/,$line);
          if(length($tem[0]) > $len)
          {
             $len=length($tem[0]);
          }
      }
      return $len;
 }

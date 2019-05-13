$num = @ARGV;
if($num != 4)
{
        die "The number of parameter is not correct!\n";
}

$input_dir = $ARGV[0];
$db_file = $ARGV[1];
$output_dir = $ARGV[2];
$scripts =  $ARGV[3];

opendir(DIR,"$input_dir") || die "Failed to open file $input_dir\n";
@files=readdir(DIR);

my($total_number)=scalar(@files);

foreach $file (@files)
{
      if($file eq "." || $file eq ".." || index($file,'network')<0)
      {
           next;
      }
      $count_number++;
      @temp = split(/\./,$file);
      $name = $temp[0];
      $networkfile = $input_dir.'/'.$file;
      print "processing $networkfile\n";
      $outputfile = $output_dir.'/'.$name.'.mitab.go_functions';
      if(-e $outputfile)
      {
          print "$outputfile has been generated, pass \n";
          next;
      }
      print("perl $scripts/A_4_prepare_go_terms_for_each_protein.pl  $networkfile $db_file $outputfile");
      $status = system("perl $scripts/A_4_prepare_go_terms_for_each_protein.pl  $networkfile $db_file $outputfile");

      if($status)
      {
        die "Failed to run <perl $scripts/A_4_prepare_go_terms_for_each_protein.pl  $networkfile $db_file $outputfile>\n"
      }

}


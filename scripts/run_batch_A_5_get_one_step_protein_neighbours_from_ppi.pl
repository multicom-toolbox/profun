$num = @ARGV;
if($num != 3)
{
        die "The number of parameter is not correct!\n";
}

$input_dir = $ARGV[0];
$output_dir = $ARGV[1];
$scripts =  $ARGV[2];

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
      $outputfile = $output_dir.'/'.$name.'.one_step_protein_neighbours';
      if(-e $outputfile)
      {
          print "$outputfile has been generated, pass \n";
          next;
      }
      print("perl $scripts/A_5_get_one_step_protein_neighbours_from_ppi.pl $networkfile $outputfile");
      $status = system("perl $scripts/A_5_get_one_step_protein_neighbours_from_ppi.pl  $networkfile $outputfile");

      if($status)
      {
        die "Failed to run <perl $scripts/A_5_get_one_step_protein_neighbours_from_ppi.pl  $networkfile $outputfile>\n"
      }

}


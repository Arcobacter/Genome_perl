#!/usr/bin/perl -w

#use Win32;

# Program BlastPTrimmer: Trims the BlastP output, eliminating matches with
# identity below a set cutoff value

@lengths = (); @taxonindex = ();
@processed1 = (); @processed2 = (); @processed3 = (); @taxa = ();
$numselected = 0; $numtaxa = 0; @labels = ''; $BLASTPname = ''; $addgenes = 1;
$currenttaxon = 0; @processed6 = ();

getinputs();
fetchtaxoninfo();
setindices();
fetchinitialtaxon();
openfiles();
while (<BLASTFILE>) {
   chomp;
   my ($Qgene,$Sgene,$ID,$VA,$VB,$VC,$QS,$QE,$SS,$SE,$Eval,$Hits) = split(/\t/);
   $VA = ''; $VB = ''; $D1 = 0; $D2 = 0; $D3 = 0; $D4 = 0; $len1 = 0; $len2 = 0;
   $Qtaxon = 0; $Flag = 0; $Staxon = 0; $gene = 0;
   ($Qtaxon,$len1) = fetchindex($Qgene);
   ($Staxon,$len2) = fetchindex($Sgene);
   $Qtype = $taxa[$Qtaxon-1]->{'Ttype'}; $Stype = $taxa[$Staxon-1]->{'Ttype'};
   $gene = parseID($Sgene,$Staxon);
   #print "$len1\t$len2\t$Qtaxon\t$Staxon\n";
   if (($Qtaxon != $currenttaxon) || eof) {
       $currenttaxon = $Qtaxon;
       sortone();
       mergesplitgenes();
       sorttwo();
       passtwo();
       printfile();
       undef(@processed1); undef(@processed2); undef(@processed3);
       undef(@processed4); undef(@processed5);
       if (eof) {last;}
       if ($taxa[$currenttaxon-1]->{'Selected'} == 1) {
          print "Processing taxon #$Qtaxon: $taxa[$Qtaxon-1]->{'Full'}\n";
       }
   }
   if (($taxa[$currenttaxon-1]->{'Selected'} == 1) && ($ID >= $simcutoff)) {
      #some gene names have an underscore character which is not good if left in
      $gene =~ s/\_//;
      #trim c and a from end of gene names, such as Cj0004c
      $gene =~ s/a$//; $gene =~ s/c$//; $gene =~ s/CDS//;
      $gene =~ s/R$//; $gene =~ s/D$//;
      #Fix the stupid W. succinogenes genome
      if ($gene eq 'pycA') {$gene = '1289';}
      $D1 = ($QE-$QS+1)/$len1; $D2 = ($QE-$QS+1)/$len2;
      $D3 = ($SE-$SS+1)/$len1; $D4 = ($SE-$SS+1)/$len2;
      $Adj = sprintf("%.1f",($Hits*(1-abs(1-$D1))*(1-abs(1-$D2))*(1-abs(1-$D3))*(1-abs(1-$D4))));
      if (($D1 < $alignlow) || ($D2 < $alignlow) || ($D3 < $alignlow) || ($D4 < $alignlow)
         || ($D1 > $alignhigh) || ($D2 > $alignhigh)|| ($D3 > $alignhigh) || ($D4 > $alignhigh)) {$Flag = 1;}
      push (@processed1,{Qgene => $Qgene, Qtaxon => $Qtaxon, Qtype => $Qtype,
                         Sgene => $Sgene, Staxon => $Staxon, Stype => $Stype, gene => $gene,
                         ID => $ID, VC => $VC, QS => $QS, QE => $QE, len1 => $len1, SS => $SS,
                         SE => $SE, len2 => $len2, Eval => $Eval, Hits => $Hits, Adj => $Adj,
                         D1 => $D1, D2 => $D2, D3 => $D3, D4 => $D4, Flag => $Flag});
   }
}
makesimmatrices();
closefiles();

sub openfiles {
   open TRIMFILE, ">Trimmed.txt";
   open BLASTFILE, "$BLASTPname";
}

sub closefiles {
   close TRIMFILE;
   close BLASTFILE;
}

sub fetchinitialtaxon {
   open BLASTFILE, "$BLASTPname";
   while (<BLASTFILE>) {
      chomp;
      my ($Qgene,$Sgene,$ID,$VA,$VB,$VC,$QS,$QE,$SS,$SE,$Eval,$Hits) = split(/\t/);
      $Qtaxon = @{$taxonindex{$Qgene}}; $Qtype = $taxa[$Qtaxon-1]->{'Ttype'};
      $currenttaxon = $Qtaxon;
      last;
   }
   print "Processing taxon #$Qtaxon: $taxa[$Qtaxon-1]->{'Full'}\n";
   close BLASTFILE;
}

sub getinputs {
    #enter output file name
    print "Enter BLASTP output file name: ";
    chomp ($BLASTPname = <STDIN>);

    #matches with %ID below $simcutoff will be excluded
    print "Enter % similarity cutoff value (default value = 30): ";
    $line = <STDIN>;
    if ($line eq "\n") {$simcutoff = 30;} else {chomp($line); $simcutoff = $line;}

    #matches with alignment lengths below alignlow or above alignhigh will be excluded;
    #for example, an alignset value of 75 removes matches whose alignment length is less
    #than 75% or more than 133% of the protein length
    print "Enter % alignment length value (default value = 75): ";
    $line = <STDIN>;
    if ($line eq "\n") {$aligncutoff = 75;} else {chomp($line); $aligncutoff = $line;}
    $alignlow = $aligncutoff/100;
    $alignhigh = 1/$alignlow;


    print "Incorporate genes into similarity matrix [Y/N]?: ";
    chomp($line = <STDIN>);
    if ($line =~ /N|n/) {$addgenes = 0; print "Gene assignment deactivated...\n";}
}

#Taxoninfo.txt contains five columns set by the user:
#column 1 is the taxon label within the gene name, e.g. CJE for CJE0001 and Ccan for Ccan_0001
#column 2 is the shorter taxon label that will be used in the printouts
#column 3 is the taxon number - this column is generally sequential from 1 to max taxa
#column 4 is the taxon type
#column 5 selects the taxa to be used in the matrices: 1 is yes, 0 is no
#column 6 is the full name/strain of each taxon

sub fetchtaxoninfo {
    print "Loading taxon information...\n";
    open TAXA, "G:/Annotations/BlastPTrimmer/Taxoninfo.txt";
    while (<TAXA>) {
       chomp;
       my ($Taxon,$Tlabel,$Tnum,$Ttype,$Selected,$Fullname) = split(/\t/);
       push (@taxa,{Taxon => $Taxon, Tlabel => $Tlabel, Tnum => $Tnum, Ttype => $Ttype, Selected => $Selected, Full => $Fullname});
       push (@labels, $Tlabel);
       $numtaxa++;
       if ($Selected == 1) {$numselected++;}
    }
    close TAXA;
}

sub setindices {
    my $sequence = "";
    my $filestart = 0;
    my $len = 0;
    my $tnum = 0;
    #my $f = 0;
    open LENGTHS, "G:/Annotations/BlastPTrimmer/Jejuniall";
    print "Loading indices...\n";
    while (<LENGTHS>) {
       chomp;
       $seg = $_;
       if ($filestart == 0) {
          $tmpname = ($_ =~ /^>(.+)/);
          $name = "$1";
          #foreach (@taxa) {if ($name =~ /$_->{'Taxon'}/) {$f = $_->{'Tnum'};}}
          #print "\tLoading taxon # $f\t$taxa[$f-1]->{'Full'}\n";
          $filestart = 1;
          next;
       }
       if ($_ =~ /^>/) {
           $len = length($sequence);
           push @{$lengths{$name}}, $len;
           foreach (@taxa) {if ($name =~ /$_->{'Taxon'}/) {$tnum = $_->{'Tnum'};}}
           push @{$taxonindex{$name}}, $tnum;
           #if ($tnum != $f) {
           #   $f = $tnum;
           #   print "\tLoading taxon # $f\t$taxa[$f-1]->{'Full'}\n";
           #}
           $tmpname = ($_ =~ /^>(.+)/);
           $name = "$1"; $tnum = 0; $sequence = "";
       } else { $sequence .= $seg; }
    }
    $len = length($sequence);
    push @{$lengths{$name}}, $len;
    foreach (@taxa) {if ($name =~ /$_->{'Taxon'}/) {$tnum = $_->{'Tnum'};}}
    push @{$taxonindex{$name}}, $tnum;
    close LENGTHS;
}

sub sortone {
    print "\tSorting BLAST file...\n";
    @processed2 = sort{$a->{Qtaxon} <=> $b->{Qtaxon} || $a->{Qgene} cmp $b->{Qgene}
                  || $a->{Staxon} <=> $b->{Staxon} || $a->{Sgene} cmp $b->{Sgene}} @processed1;
}

sub sorttwo {
    print "\tSorting BLAST file...\n";
    @processed4 = sort{$a->{Qtaxon} <=> $b->{Qtaxon} || $a->{Qgene} cmp $b->{Qgene}
                  || $a->{Staxon} <=> $b->{Staxon} || $b->{Adj} <=> $a->{Adj}} @processed3;
}

sub passtwo {
    print "\tProcessing BLAST file...Pass 2\n";
    $start = 1; $setbug1 = ''; $setbug2 = 0;
    foreach (@processed4) {
       if ($start == 1) {
          if (($_->{'Flag'} == 0) || ($_->{'Flag'} == 2) || ($_->{'Flag'} == 4)) {
              push @processed5, $_;
          }
          $start = 0;
          $setbug1 = $_->{'Qgene'};
          $setbug2 = $_->{'Staxon'};
       }
       if ($start == 0) {
          if ($_->{'Qgene'} ne $setbug1) {
              if (($_->{'Flag'} == 0) || ($_->{'Flag'} == 2) || ($_->{'Flag'} == 4)) {
                  push @processed5, $_;
              }
              $setbug1 = $_->{'Qgene'};
              $setbug2 = $_->{'Staxon'};
          }
          if ($_->{'Qgene'} eq $setbug1) {
             if ($_->{'Staxon'} != $setbug2) {
                 $setbug2 = $_->{'Staxon'};
                 if (($_->{'Flag'} == 0) || ($_->{'Flag'} == 2) || ($_->{'Flag'} == 4)) {
                     push @processed5, $_;
                 }
             }
          }
       }
    }
}

sub parseID {
    my($genenum,$bugnum) = @_;
    foreach (@taxa) {
       if ($bugnum == $_->{'Tnum'})  {$genenum =~ s/$_->{'Taxon'}//;}
    }
    return ($genenum);
}

sub printrecord {
    my($field) = @_;
    my($num) = 0;
    $num = $field->{'Qtaxon'};
    @temparray = ();
    push (@temparray, (
          $field->{'Qgene'},  $field->{'Qtaxon'}, $field->{'Qtype'},
          $field->{'Sgene'},  $field->{'gene'},   $field->{'Staxon'},
          $field->{'Stype'},  $field->{'ID'},
          $field->{'VC'},     $field->{'QS'},     $field->{'QE'},
          $field->{'len1'},   $field->{'SS'},     $field->{'SE'},
          $field->{'len2'},   $field->{'Eval'},   $field->{'Hits'},
          $field->{'Adj'},    $field->{'D1'},     $field->{'D2'},
          $field->{'D3'},     $field->{'D4'},     $field->{'Flag'}));
    $tabbedline = join("\t",@temparray);
    print TRIMFILE "$tabbedline\n";
}

sub fetchindex  {
   my ($g) = @_;
   @Len = @{$lengths{$g}};
   @Tax = @{$taxonindex{$g}};
   return (@Tax,@Len);
}

sub printfile {
    print "\tPrinting BLAST file...\n";
    foreach (@processed5) {
       printrecord($_);
       push (@processed6,{Qgene => $_->{'Qgene'},
                          Qtaxon => $_->{'Qtaxon'},
                          Sgene => $_->{'Sgene'},
                          Staxon => $_->{'Staxon'},
                          ID => $_->{'ID'},});
    }
}

sub mergesplitgenes {
   print "\tMerging split genes...\n";
   $slot1 = 0; $slot2 = 0; $start = 1; $split = 0;
   foreach (@processed2) {
      if ($start == 1) {$start = 0; $slot1 = $_; next;}
      if ($start == 0) {
         if ($slot2 == 0) {$slot2 = $_; next;}
         if (($slot2->{'gene'} == $slot1->{'gene'} + 1)
            && ($slot1->{'Staxon'} eq $slot2->{'Staxon'})
            && ($slot1->{'D1'} + $slot2->{'D1'} >= $alignlow)
            && ($slot1->{'D1'} + $slot2->{'D1'} <= $alignhigh)
            && ($slot1->{'D3'} + $slot2->{'D3'} >= $alignlow)
            && ($slot1->{'D3'} + $slot2->{'D3'} <= $alignhigh)
            && (($slot2->{'QS'} >= $slot1->{'QE'}) || ($slot1->{'QS'} > $slot2->{'QE'}))) {
              $slot2->{'Sgene'} = $slot1->{'Sgene'} . "\/" . $slot2->{'gene'};

              $slot2->{'ID'} = ($slot2->{'ID'} + $slot1->{'ID'})/2;
              $slot2->{'len2'} = $slot1->{'len2'} + $slot2->{'len2'};

              if ($slot1->{'QS'} < $slot2->{'QS'}) {$slot2->{'QS'} = $slot1->{'QS'};}
              else {$slot2->{'QE'} = $slot1->{'QE'}; $slot2->{'SS'} = $slot1->{'SS'};}

              $slot2->{'SE'} = $slot1->{'SE'} + $slot2->{'SE'};
              $slot2->{'Hits'} = $slot1->{'Hits'} + $slot2->{'Hits'};
              if ($slot1->{'Eval'} < $slot2->{'Eval'}) {$slot2->{'Eval'} = $slot1->{'Eval'};}
              $slot2->{'D1'} = ($slot2->{'QE'}-$slot2->{'QS'}+1)/$slot2->{'len1'};
              $slot2->{'D2'} = ($slot2->{'QE'}-$slot2->{'QS'}+1)/$slot2->{'len2'};
              $slot2->{'D3'} = ($slot2->{'SE'}-$slot2->{'SS'}+1)/$slot2->{'len1'};
              $slot2->{'D4'} = ($slot2->{'SE'}-$slot2->{'SS'}+1)/$slot2->{'len2'};
              $slot2->{'Adj'} = $slot2->{'Hits'}*(1-abs(1-$slot2->{'D1'}))*(1-abs(1-$slot2->{'D2'}))*(1-abs(1-$slot2->{'D3'}))*(1-abs(1-$slot2->{'D4'}));
              $slot1->{'Flag'} = 99;
              $slot2->{'Flag'} = 3;
              if (($slot2->{'D1'} > $alignlow) && ($slot2->{'D2'} > $alignlow) && ($slot2->{'D3'} > $alignlow) && ($slot2->{'D4'} > $alignlow)
                 && ($slot2->{'D1'} < $alignhigh) && ($slot2->{'D2'} < $alignhigh) && ($slot2->{'D3'} < $alignhigh) && ($slot2->{'D4'} < $alignhigh)) {$slot2->{'Flag'} = 2}
              $slot1 = $slot2;
              $slot2 = $_;
              $split++;
              next;
         }
         if (($slot2->{'gene'} == $slot1->{'gene'})
            && ($slot1->{'Sgene'} eq $slot2->{'Sgene'})
            && ($slot2->{'len1'} / $slot2->{'len2'} <= $alignhigh)
            && ($slot2->{'len2'} / $slot2->{'len1'} <= $alignhigh)
            && ($slot1->{'Flag'} == 1)
            && ($slot2->{'Flag'} == 1)) {
              $slot2->{'ID'} = ($slot2->{'ID'} + $slot1->{'ID'})/2;
              if ($slot1->{'QS'} < $slot2->{'QS'}) {$slot2->{'QS'} = $slot1->{'QS'};} else {$slot2->{'QE'} = $slot1->{'QE'};}
              if ($slot1->{'SS'} < $slot2->{'SS'}) {$slot2->{'SS'} = $slot1->{'SS'};} else {$slot2->{'SE'} = $slot1->{'SE'};}
              $slot2->{'Hits'} = $slot1->{'Hits'} + $slot2->{'Hits'};
              if ($slot1->{'Eval'} < $slot2->{'Eval'}) {$slot2->{'Eval'} = $slot1->{'Eval'};}
              $slot2->{'D1'} += $slot1->{'D1'};
              $slot2->{'D2'} += $slot1->{'D2'};
              $slot2->{'D3'} += $slot1->{'D3'};
              $slot2->{'D4'} += $slot1->{'D4'};
              $slot2->{'Adj'} = $slot2->{'Hits'}*(1-abs(1-$slot2->{'D1'}))*(1-abs(1-$slot2->{'D2'}))*(1-abs(1-$slot2->{'D3'}))*(1-abs(1-$slot2->{'D4'}));
              $slot1->{'Flag'} = 99;
              $slot2->{'Flag'} = 3;
              if (($slot2->{'D1'} > $alignlow) && ($slot2->{'D2'} > $alignlow) && ($slot2->{'D3'} > $alignlow) && ($slot2->{'D4'} > $alignlow)
                 && ($slot2->{'D1'} < $alignhigh) && ($slot2->{'D2'} < $alignhigh) && ($slot2->{'D3'} < $alignhigh) && ($slot2->{'D4'} < $alignhigh)) {$slot2->{'Flag'} = 4}
              $slot1 = $slot2;
              $slot2 = $_;
              $split++;
              next;
         }
         else {
            if ($slot1->{'Flag'} != 99) {push @processed3, $slot1;}
            $slot1 = $slot2;
            $slot2 = $_;
         }
      }
   }
   push @processed3, $slot2;
   print "\t$split genes processed\n";
}

sub makesimmatrices {
   my @identities = '';
   my @Sgenes = '';
   my @coregenes = '';
   my @specificgenes = '';
   my @homologs = '';
   my @novelgenes = ();
   my $genecounter = 0;
   my $start = 1;
   my $startgene = '';
   my $novel = 1;
   print "Constructing similarity matrices...\n";
   open MATRIXFILE, ">Genesimmatrix.txt";
   open ORTHOLOGFILE, ">Orthologlist.txt";
   print MATRIXFILE "Qgene";
   if ($addgenes == 0) {
      for ($i = 1; $i <= $numtaxa; $i++) {print MATRIXFILE "\t$labels[$i]";}
   }
   if ($addgenes == 1) {
      for ($i = 1; $i <= $numtaxa; $i++) {print MATRIXFILE "\t\t$labels[$i]";}
   }
   print MATRIXFILE "\n";
   for ($i = 1; $i <= $numtaxa; $i++) {
      #pairwise identity scores and subject genes
      $identities[$i] = 0; $Sgenes[$i] = '';
      #genes conserved across all selected taxa
      $coregenes[$i] = 0;
      #genes specific to selected taxa
      $specificgenes[$i] = 0;
      #computes the size of the gene pool for selected taxa
      $novelgenes[$i] = 0;
      #compiles identities for Taxonsimmatrix
      for ($j = 1; $j <= $numtaxa; $j++) {$homologs[$i][$j] = 0;}
   }
   foreach (@processed6) {
      if ($start == 1) {
         $startgene = $_->{'Qgene'};
         $identities[$_->{'Staxon'}] = $_->{'ID'};
         $Sgenes[$_->{'Staxon'}] = $_->{'Sgene'};
         $start = 0;
         #first gene - gene pool size is 1
         $novelgenes[$_->{'Qtaxon'}]++;
      }
      if ($start == 0) {
         #same gene - add ID score (and gene) to array
         if ($startgene eq $_->{'Qgene'}) {
            $identities[$_->{'Staxon'}] = $_->{'ID'};
            $Sgenes[$_->{'Staxon'}] = $_->{'Sgene'};
         }
         #next gene in list
         if ($startgene ne $_->{'Qgene'}) {
             #print ID scores for previous gene
             print MATRIXFILE "$startgene\t";
             print ORTHOLOGFILE "$startgene\t";
             for ($i = 1; $i <= $numtaxa; $i++) {
                if ($addgenes == 0) {
                   if ($identities[$i] == 0) {print MATRIXFILE "---\t";}
                   else {
                      print MATRIXFILE "$identities[$i]\t";
                      print ORTHOLOGFILE "$Sgenes[$i]\t";
                   }
                }
                if ($addgenes == 1) {
                   if ($identities[$i] == 0) {print MATRIXFILE "---\t---\t";}
                   else {
                      print MATRIXFILE "$Sgenes[$i]\t$identities[$i]\t";
                      print ORTHOLOGFILE "$Sgenes[$i]\t";
                   }
                }
             }

             #within selected set how many genes showed a positive match (%ID>0)
             for ($i = 1; $i <= $numtaxa; $i++) {
                if (($identities[$i] > 0) && ($taxa[$i-1]->{'Selected'} == 1)) {$genecounter++;}
             }
             #if only one gene matches (self match) then increment specific gene counter
             if ($genecounter == 1) {$specificgenes[$_->{'Qtaxon'}]++;}
             #if all genes in selected set match, gene is a core gene
             if ($genecounter == $numselected) {
                $coregenes[$_->{'Qtaxon'}]++;
                for ($i = 1; $i <= $numtaxa; $i++) {
                   $homologs[$_->{'Qtaxon'}][$i] += $identities[$i];
                }
             }

             if ($_->{'Qtaxon'} == 1) {$novelgenes[1]++;}
             else {
                if ($_->{'Qtaxon'} <= $numtaxa) {
                   for ($i = 1; $i < $_->{'Qtaxon'}; $i++) {
                      if (($identities[$i] > 0) && ($taxa[$i-1]->{'Selected'} == 1)) {
                         $novel = 0;
                      }
                   }
                   if ($novel == 1) {$novelgenes[$_->{'Qtaxon'}]++;}
                   $novel = 1;
                }
             }

             print MATRIXFILE "\n";
             print ORTHOLOGFILE "\n";
             $startgene = $_->{'Qgene'};
             for ($i = 1; $i <= $numtaxa; $i++) {$identities[$i] = 0; $Sgenes[$i] = 0;}
             $identities[$_->{'Staxon'}] = $_->{'ID'};
             $Sgenes[$_->{'Staxon'}] = $_->{'Sgene'};
             $genecounter = 0;
         }
      }
   }
   print MATRIXFILE "$startgene\t";
   for ($i = 1; $i <= $numtaxa; $i++) {
      if ($addgenes == 0) {
         if ($identities[$i] == 0) {print MATRIXFILE "---\t";}
         else {
            print MATRIXFILE "$identities[$i]\t";
            print ORTHOLOGFILE "$Sgenes[$i]\t";
         }
      }
      if ($addgenes == 1) {
         if ($identities[$i] == 0) {print MATRIXFILE "---\t---\t";}
         else {
            print MATRIXFILE "$Sgenes[$i]\t$identities[$i]\t";
            print ORTHOLOGFILE "$Sgenes[$i]\t";
         }
      }
   }
   close MATRIXFILE;
   close ORTHOLOGFILE;

   open TYPES, ">Genedata.txt";
   print TYPES "Taxon\tNovel\tSpecific\tCore\n";
   for ($i = 1; $i <= $numtaxa; $i++) {
      if ($taxa[$i-1]->{'Selected'} == 1) {print TYPES "$labels[$i]\t";}
      if ($taxa[$i-1]->{'Selected'} == 1) {print TYPES "$novelgenes[$i]\t";}
      if ($taxa[$i-1]->{'Selected'} == 1) {print TYPES "$specificgenes[$i]\t";}
      if ($taxa[$i-1]->{'Selected'} == 1) {print TYPES "$coregenes[$i]\n";}
   }
   close TYPES;

   # create pairwise taxon vs. taxon aa identity grid
   open HOMOLOGFILE, ">Taxonsimmatrix.txt";
   for ($i = 1; $i <= $numtaxa; $i++) {
     if ($taxa[$i-1]->{'Selected'} == 1) {print HOMOLOGFILE "\t$labels[$i]";}
   }
   print HOMOLOGFILE "\n";
   for ($i = 1; $i <= $numtaxa; $i++) {
      if ($taxa[$i-1]->{'Selected'} == 1) {
         print HOMOLOGFILE "$labels[$i]\t";
         for ($j = 1; $j <= $numtaxa; $j++) {
            if ($taxa[$j-1]->{'Selected'} == 1) {
               if ($coregenes[$i] > 0) {
                  # round aa identity values to nearest integer
                  $ID = int(($homologs[$i][$j] / $coregenes[$i]) + 0.5);
               }
               else {$ID = 0;}
               print HOMOLOGFILE "$ID\t";
            }
         }
      print HOMOLOGFILE "\n";
      }
   }
   close HOMOLOGFILE;
}
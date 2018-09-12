#!/usr/bin/perl -w
use Win32;

$sequence = ''; # sequence as one long string
$seqlen = 0; # length of sequence
@contigs = (); # array of hash refs with name/length/sequence keys
@sffreads = (); # array of hash refs with name/length/sequence keys
$filestart = 0;
@reads = (); # array of hash refs with name/position/sequence keys
$contigID = '';

# initialization parameters

$contigfile = "q:/454 Genomes/C_seal_17262/454AllContigs.fna";
$readfile = "q:/454 Genomes/Sff_fasta/17262RL8.fas";
$outputfile = "q:/454 Genomes/Flanking/17262Flanking.txt";
$endgrab = 15;    #nt window at contig ends to fish out overlapping reads
$terse = 0; #($terse=0) = show only reads with defined matches, ($terse=1) = shows all reads

# grab contigs and load sequences into hash array
loadcontigs(); $filestart = 0;

# grab SFF reads and load sequences into hash array
loadsfffile(); $filestart = 0;

@contigs2 = @contigs;
open(OUTPUT,">$outputfile");
$Lcontigend = ''; $Rcontigend = '';
foreach (@contigs2) {
   #grab 15 base terminal sequences of each contig, left and right ends
   $Lcontigend = substr($_->{'data'},0,$endgrab);
   $Rcontigend = substr($_->{'data'},$_->{'length'}-$endgrab,$endgrab);
   #convert sequences to upper case
   $Lcontigend =~ tr/a-z/A-Z/;
   $Rcontigend =~ tr/a-z/A-Z/;
   $contigID = $_->{'ID'};
   print "Processing $contigID left\n";
   #find left-adjoining contigs
   grabreads($Lcontigend,"L");
   findcontigs("L");
   @reads = ();
   print "Processing $contigID right\n";
   #find right-adjoining contigs
   grabreads($Rcontigend,"R");
   findcontigs("R");
   @reads = ();

   print OUTPUT "\n";
}
close(OUTPUT);

sub loadcontigs {
   print "Loading contig sequences...\n";
   open(FASTA,"$contigfile") or die "couldn't open sequence_file: $!\n";
   while (<FASTA>) {
      chomp;
      $seg = $_;
      if ($filestart == 0) {
         $_ =~ /^>(.+)  len/;
         $seqID = $1;
         $filestart = 1;
         next;
      }
      if ($_ =~ /^>(.+)  len/) {
         $seqlen = length $sequence;
         $sequence =~ tr/a-z/A-Z/;
          push(@contigs, {ID => $seqID, length => $seqlen, data => $sequence});
         $seqID = $1; $sequence = '';
      } else { $sequence .= $seg; }
   }
   $seqlen = length $sequence;
   push(@contigs, {ID => $seqID, length => $seqlen, data => $sequence});
   close(FASTA);
}

sub loadsfffile {
   print "Loading sff reads into memory...\n";
   $sffnum = 0;
   open(SFFFILE,"$readfile") or die "couldn't open sequence_file: $!\n";
   while (<SFFFILE>) {
      chomp;
      $seg = $_;
      if ($filestart == 0) {
         $_ =~ /^>(.+) len/;
         $sffID = $1;
         $filestart = 1;
         next;
      }
      if ($_ =~ /^>(.+) len/) {
         $seqlen = length $sequence;
         $sffnum++;
         push(@sffreads, {ID => $sffID, length => $seqlen, data => $sequence});
         $sffID = $1; $sequence = '';
      } else { $sequence .= $seg; }
   }
   $sffnum++;
   push(@sffreads, {ID => $sffID, length => $seqlen, data => $sequence});
   close(SFFFILE);
   print "$sffnum reads loaded into memory\n";
}

sub grabreads {
   my ($contigend,$orient) = @_;
   $readOK = 0; $loc = 0;
   foreach $s (@sffreads) {
      if (($s->{'data'} =~ /$contigend/ig) && ($s->{'length'} >= 250)) {
         $readOK = 1; $loc = pos($s->{'data'})-14;
      }
      if (($loc <= 15) && ($orient eq "L")) {$readOK = 0;}
      if (($loc >= $s->{'length'} - 30) && ($orient eq "R")) {$readOK = 0;}
      if ($readOK == 1) {push(@reads, {ID => $s->{'ID'}, loc => $loc, length => $s->{'length'}, data => $s->{'data'}});}
      $readOK = 0; $loc = 0;

      $compread = '';
      $compread = reverse($s->{'data'});
      $compread =~ tr/ACGTacgt/TGCAtgca/;
      if (($compread =~ /$contigend/ig) && ($s->{'length'} >= 250)) {
         $readOK = 1; $loc = pos($compread)-14;
      }
      if (($loc <= 15) && ($orient eq "L")) {$readOK = 0;}
      if (($loc >= $s->{'length'} - 30) && ($orient eq "R")) {$readOK = 0;}
      if ($readOK == 1) {
          $IDc = $s->{'ID'} . "c";
          push(@reads, {ID => $IDc, loc => $loc, length => $s->{'length'}, data => $compread});
      }
      $readOK = 0; $loc = 0;
   }
}

sub findcontigs {
   my ($orient) = @_;
   $query = '';
   if ($orient eq "L") {
      foreach $i (@reads) {
         $ended = 0;
         $start = $i->{'loc'}-16; $multiple = 0; $numcontig = '';
         until (($start == 0) || ($ended == 1)) {
            $query = substr($i->{'data'},$start,15);
            foreach $j (@contigs) {
               $revquery = reverse $query;
               $revquery =~ tr/AGCT/TCGA/;
               if ($j->{'data'} =~ /$query$/i) {
                  if ($start - $j->{'length'} < 0) {$ended = 1;}
                  if ($multiple == 1) {$numcontig .= ", ";}
                  $numcontig .= $j->{'ID'};
                  $multiple = 1;
               }
               if ($j->{'data'} =~ /^$revquery/i) {
                  if ($start - $j->{'length'} < 0) {$ended = 1;}
                  if ($multiple == 1) {$numcontig .= ", ";}
                  $numcontig .= $j->{'ID'};
                  $numcontig .= "c";
                  $multiple = 1;
               }
            }
            $start--;
         }
         if ($terse == 0) {
           if ($numcontig ne '') {
              print OUTPUT "$contigID\t$i->{'ID'}\t$i->{'loc'}\t$i->{'length'}\t$orient\t$numcontig\n";
           }
         }
         if ($terse == 1) {
            print OUTPUT "$contigID\t$i->{'ID'}\t$i->{'loc'}\t$i->{'length'}\t$orient\t$numcontig\n";
         }
         $query = '';
      }
   }
   if ($orient eq "R") {
      foreach $i (@reads) {
         $ended = 0;
         $start = $i->{'loc'}+14; $multiple = 0; $numcontig = '';
         until (($start >= $i->{'length'}-15) || ($ended ==1)) {
            $query = substr($i->{'data'},$start,15);
            foreach $j (@contigs) {
               $revquery = reverse $query;
               $revquery =~ tr/AGCT/TCGA/;
               if ($j->{'data'} =~ /^$query/i) {
                  if ($start + $j->{'length'} > $i->{'length'}) {$ended = 1;}
                  if ($multiple == 1) {$numcontig .= ", ";}
                  $numcontig .= $j->{'ID'};
                  $multiple = 1;
               }
               if ($j->{'data'} =~ /$revquery$/i) {
                  if ($start + $j->{'length'} > $i->{'length'}) {$ended = 1;}
                  if ($multiple == 1) {$numcontig .= ", ";}
                  $numcontig .= $j->{'ID'};
                  $numcontig .= "c";
                  $multiple = 1;
               }
            }
            $start++;
          }
         if ($terse == 0) {
           if ($numcontig ne '') {
              print OUTPUT "$contigID\t$i->{'ID'}\t$i->{'loc'}\t$i->{'length'}\t$orient\t$numcontig\n";
           }
         }
         if ($terse == 1) {
            print OUTPUT "$contigID\t$i->{'ID'}\t$i->{'loc'}\t$i->{'length'}\t$orient\t$numcontig\n";
         }
         $query = '';
      }
   }
}
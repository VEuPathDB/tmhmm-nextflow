#!/usr/bin/perl

use strict;

use Getopt::Long;

use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;


my ($file, $outputFile);

&GetOptions(
    "inputFile=s"       => \$file,
    "outputFile=s"       => \$outputFile,
    );


open(IN, "<$file") or die "Couldn't open file '$file': $!\n";

open(OUT, ">$outputFile") or die "Couldn't open file $outputFile for writing: $!";
my $source = "veupathdb";


while (<IN>) {
    chomp;

    my ($attVals, $locs) =  &parseAttVals($_);
    next if (!$attVals);

    my $sourceId = $attVals->{protein_id};

    my $feat = Bio::SeqFeature::Generic->new(
        -start        => 1,
        -end          => $attVals->{'length'},
        -primary      => "tmhmm2.0", # -primary_tag is a synonym
        -seq_id       => $sourceId,
        -source_tag   => $source,
        -score        => ".",
        -frame        => ".",
        -tag          => {ID => "${sourceId}_tmhmm",
                          ExpectedAA => $attVals->{'expected_aa'},
                          First60 =>   $attVals->{'first_60'} ,
                          PredictedHelices => $attVals->{'predicted_helices'}
                        } );


    $feat->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
    print OUT $feat->gff_string, "\n" ;

    my $count = 1;

    foreach my $loc (@$locs) {
        my $type = $loc->[0];
        my $start = $loc->[1];
        my $end = $loc->[2];

        if ($type && $start && $end) {
            my $subfeat = Bio::SeqFeature::Generic->new(
                -start        => $start,
                -end          => $end,
                -primary      => $type, # -primary_tag is a synonym
                -seq_id       => $sourceId,
                -source_tag   => $source,
                -score        => ".",
                -frame        => ".",
                -tag          => {ID => "${sourceId}_tmhmm_${count}",
                                  Parent => "${sourceId}_tmhmm"

                } );


            $subfeat->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
            print OUT $subfeat->gff_string, "\n" ;
        }
        $count++;
    }
}


close OUT;
close IN;

sub parseAttVals{
  my ($line_in) = @_;

  $line_in =~ s/\w+\=//g;

  my %attVals;

  my @Record = split(/\t/,$line_in);

  return if ($Record[4] == 0);

  $attVals{'protein_id'} = $Record[0];

  $attVals{'length'} = $Record[1];

  $attVals{'expected_aa'} = $Record[2];

  $attVals{'first_60'} = $Record[3];

  $attVals{'predicted_helices'} = $Record[4];

  $attVals{'is_predicted'} = 1;

  my $topology = $Record[5];

  my ($locs) = &parseTopology($topology, $Record[1]);

  die "Can't get the TM protein type or locations for $attVals{'aa_sequence_id'}" if (scalar @$locs == 0);

  return (\%attVals,$locs);
}

sub parseTopology {
  my ($topology, $length) = @_;

  my @rv;

  $topology =~ s/\s//g;

  my @locs = split(/([io])/,$topology);

  # @locs will contain alternating start-end and delimiters (i|o)
  for (my $i = 0; $i <= $#locs; $i++) {

      my $loc = $locs[$i];

      my $prevLoc = $locs[$i-1];
      my $nextLoc = $locs[$i+1];

      next unless($loc);

      if($loc eq 'i' || $loc eq 'o') {
          my $start = $prevLoc ? &getEndLoc($prevLoc) + 1 : 1;
          my $end = $nextLoc ? &getStartLoc($nextLoc) - 1 : $length;

          my $type = $loc eq 'i' ? "inside" : "outside";

          push(@rv, [$type, $start, $end]);
      }

      else {
          my $loc = &parseLoc($loc);
          push @rv, ["TMhelix", $loc->[0], $loc->[1]];
      }
  }

  return \@rv;

}


sub parseLoc {
    my ($loc) = @_;

    my ($start, $end) = split(/\-/, $loc);

    return [$start, $end];

}

sub getStartLoc {
    my ($loc) = @_;

    return &parseLoc($loc)->[0];

}
sub getEndLoc {
    my ($loc) = @_;
    return &parseLoc($loc)->[1];

}


1;

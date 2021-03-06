=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

  PrimateAI

=head1 SYNOPSIS

  mv PrimateAI.pm ~/.vep/Plugins

  ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2.tsv.bgz
  ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2_hg38.tsv.bgz

=head1 DESCRIPTION

  A VEP plugin designed to retrieve clinical impact scores of variants, as described in https://www.nature.com/articles/s41588-018-0167-z.
  In brief, common missense mutations in non-human primate species are usually benign in humans. They can therefore be eliminated from screens for pathogenic mutations.

  Files containing predicted pathogenicity scores can be downloaded from https://basespace.illumina.com/s/yYGFdGih1rXL (a free BaseSpace account may be required):
      PrimateAI_scores_v0.2.tsv.gz (for GRCh37/hg19)
      PrimateAI_scores_v0.2_hg38.tsv.gz (for GRCh38/hg38)

  Before running the plugin for the first time, the following steps must be taken to format the downloaded files:
  1.  Unzip the score files
  2.  Add '#' in front of the column description line
  3.  Remove any empty lines.
  4.  Remove 'chr' prefix from data lines
  5.  Sort the file by chromosome and position
  6.  Compress the file in .bgz format
  7.  Create tabix index (requires tabix to be installed).
  Command line example:
      $ gunzip PrimateAI_scores_v0.2.tsv.gz | sed '12s/.*/#&/' | sed s/^chr// | sed '/^$/d | | awk 'NR<13{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' | bgzip > PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz
      $ tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz

=cut

package PrimateAI;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $file = $self->params->[0];
  die("ERROR: PrimateAI scores file $file not found\n") unless $file && -e $file;

  $self->add_file($file);

  #Check that the file matches the assembly
  if ($file =~ /PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz/){
    die "The PrimateAI scores file contains GRCh37 coordinates, but the assembly used by VEP is not GRCh37\n" unless $self->{config}->{assembly} eq "GRCh37";
  }

  if ($file =~ /PrimateAI_scores_v0.2i_hg38_sorted.tsv.bgz/){
    die "The PrimateAI scores file contains GRCh38 coordinates, but the assembly used by VEP is not GRCh38\n" unless $self->{config}->{assembly} eq "GRCh38";
  }

  $self->expand_left(0);
  $self->expand_right(0);

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  return {
    PrimateAI => "PrimateAI score for variants"
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $allele = $tva->variation_feature_seq;
  #my $alt_aa = $tva->peptide;

  return {} unless $allele =~ /^[ACGT]$/;

  #Get the start and end coordinates, and ensure they are the right way round (i.e. start < end).
  my $vf_start = $vf->{start};
  my $vf_end = $vf->{end};
  ($vf_start, $vf_end) = ($vf_end, $vf_start) if $vf_start > $vf_end;

  #Check the strands and complement the allele if necessary.
  if ($vf->{strand} <0){
    reverse_comp(\$allele);
  }

  #Compare the postion, allele and alt amino acid
    my ($res) = grep {
    $_->{start} eq $vf_start &&
    $_->{end} eq $vf_end &&
    $_->{alt} eq $allele
    #$_->{alt_aa} eq $alt_aa
  } @{$self->get_data($vf->{chr}, $vf_start, $vf_end)};

  #Return data if matched
  return $res ? $res->{result} : {};
}

sub parse_data {
  my ($self, $line) = @_;

  #Columns in order of the input file.
  my ($chr, $pos, $ref, $alt, $refAA, $altAA, $strand, $trinuc, $gene, $cov, $score) = split(/\t/, $line);

  return {
    chr => $chr,
    alt => $alt,
    start => $pos,
    end => $pos,
    strand => $strand,
    alt_aa => $altAA,
    result => {
      #Score to be returned
      PrimateAI   => $score

      #Are any other features useful to return?
    }
  };
}

sub get_start {
  return $_[1]->{'start'};
}

sub get_end {
  return $_[1]->{'end'};
}

1;

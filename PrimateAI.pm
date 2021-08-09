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

  ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2.tsv.gz
  ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2_hg38.tsv.gz

=head1 DESCRIPTION

  A VEP plugin designed to retrieve clinical impact scores of variants, as described in https://www.nature.com/articles/s41588-018-0167-z.
  Briefly, common missense mutations in non-human primate species are usually benign in humans, so can be eliminated from screens for pathogenic mutations.

  Files containing predicted pathogenicity scores can be downloaded from https://basespace.illumina.com/s/yYGFdGih1rXL (a free BaseSpace account may be required):
      PrimateAI_scores_v0.2.tsv.gz (for GRCh37/hg19)
      PrimateAI_scores_v0.2_hg38.tsv.gz (for GRCh38/hg38)

  Before running the plugin for the first time, the following steps must be taken to format the downloaded files:
  1.  Unzip the scores files
  2.  Add '#' in front of the column description line
  3.  Remove 'chr' prefix from data lines
  4.  Compress the file again
  5.  Create tabix index (requires tabix to be installed).
  Command line example:
      $ gunzip PrimateAI_scores_v0.2.tsv.gz | sed '12s/.*/#&/' | sed s/^chr// | bgzip > PrimateAI_scores_v0.2.tsv.bgz'
      $ tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2.tsv.bgz'

=cut

Package PrimateAI;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $file = $self->params->[0];
  $self->add_file($file);

  die("ERROR: PrimateAI scores file $file not found\n") unless $file && -e $file;

  #Check that the file matches the assembly
  if ($self->{config}->{assembly} eq "GRCh37"){
    if ($file =~ /hg38/){die "The assembly used is GRCh37, but the Primate AI file contains GRCh38/hg38 coordinates.\n";}
  } elsif ($self->{config}->{assembly} eq "GRCh38"){
      if ($file !~ /hg38/){die "The assembly used is GRCh38, but the Primate AI file contains GRCh37 coordinates.\n";}
  }

  $self->expand_left(0);
  $self->expand_right(0);

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub variation_feature_types {
  return ['VariationFeature'];
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
  my $alt_aa = $tva->peptide;

  #Check the strands and complement the allele if necessary.
  if ($vf->{strand} <0){
    #1 indicates positive strand
    if ($strand == 1){reverse_comp(\$allele);}
  } elsif{
      #0 indicates negative strand
      if ($strand == 0){reverse_comp(\$allele);}
  }

  return {} unless $allele =~ /^[ACGT]$/;

  #Get the start and end coordinates, and ensure they are the right way round (i.e. start < end).
  my $vf_start = $vf->{start};
  my $vf_end = $vf->{end};
  ($vf_start, $vf_end) = ($vf_end, $vf_start) if $vf_start > $vf_end;

#Compare the postion, allele and alt amino acid
  my ($res) = grep {
    $_->{pos} eq $vf_start &&
    $_->{pos} eq $vf_end &&
    $_->{alt} eq $allele &&
    $_->{altAA} eq $alt_aa
  } @{$self->get_data($vf->{chr}, $vf_start, $vf_end)};

  #Return data if matched
  return $res ? $res->{result} : {};
}

sub parse_data {
  my ($self, $line) = @_;

  #Columns in order of the input file.
  my ($chr, $pos, $ref, $alt, $refAA, $altAA, $strand, $trinuc, $gene, $cov, $score) = split(/\t/, $line);

  return {
    ref => $chr,
    alt => $alt,
    start => $pos,
    end => $pos,
    result => {
      #Score to be returned
      PrimateAI   => $score,

      #Are any other features useful to return?
    }
  };
}

sub get_start {
  return $_[1]->{'pos'};
}

sub get_end {
  return $_[1]->{'pos'};
}

1;

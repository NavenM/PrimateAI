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

=head1 DESCRIPTION
 
  A VEP plugin designed to retrieve clinical impact scores of variants, as described in https://www.nature.com/articles/s41588-018-0167-z.

=cut

Package PrimateAI

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  #Check that tabix is present and correct.
  die "Warning: tabix not found.\n"; unless `which tabix 2>$1` =~ /tabix$/;
  
  my $file = $self->params->[0];

  die("ERROR: PrimateAI scores file $file not found\n") unless $file && -e $file;

  #$self->expand_left(0);
  #$self->expand_right(0);

  open FILE, "$file" or die "File ($file) could not be opened. Exiting.\n";
  while (<FILE>){
    chomp;
    
    next unless $file =~ /^chr\t/;

    $self->{headers} = [split];
  }

  close FILE;
  
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
  
  return {} unless $allele =~ /^[ACGT]+$/;

  my @data =  @{$self->get_data($vf->{chr}, $vf->{start}, $vf->{end})};

  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;

  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;

  my ($res) = grep {
    $_->{chr} eq $vf->{chr} &&
    $_->{pos} eq $vf->{start} &&
    $_->{pos} eq $vf->{end} &&
    $_->{alt} eq $allele
  } @{$self->get_data($vf->{chr}, $vf->{start}, $vf->{end})};

  return {} unless(@data);
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $pos, $ref, $alt, $refAA, $altAA, $strand, $trinuc, $gene, $cov, $score) = split /\t/, $line;

  return {
    ref => $chr,
    alt => $alt,
    start => $pos,
    end => $pos,
    result => {
      PrimateAI   => $score,
    }
  };
}

sub get_start {  
  return $_[1]->{'pos'};
}

sub get_end {
  return $_[1]->{'pos'};
}






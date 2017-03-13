#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use IO::File;
use constant DEBUG=>0;

=head1 vcf_vetting.pl 
 
 This script is based on earlier work on merging mutect/strelka files
 and essentially performs the following tasks:

 * removes non-canonical contigs

 * adds GT field (dot or calculated based on SGT, if available)

 * removes tool-specific header lines

 this script will take the strelka or mutect vcf files and clean the headers, also removing non-canonical contigs
 will impute GT field if it is absent

=cut

my $USAGE = "vcf_vetting.pl --input [input vcf file (may be bgzipped)] --infiles [optional comma-sep list of all used inputs] --index [optional index to use for disambiguation, needed only when infiles are passed] --output [output file] \n";
my %opts = (canonical=>1);  ### only show canonical

GetOptions('input=s'              => \$opts{input},
           'source=s'             => \$opts{source},
           'infiles|inlist|all=s' => \$opts{infiles},
           'index|suffix=i'       => \$opts{index},
           'output|out=s'         => \$opts{output});
my $input  = $opts{input};
my $output = $opts{output};
if (! $input || !-e $input || !$output || ($opts{infiles} && !$opts{index})) { die $USAGE; }

### Check if we have an input

die "Need some input to work with" if ! -e $input;

### Check if we have infiles and if yes, there are two of them
my $formats = {};

if ($opts{infiles}) {
 my @ins = split(",", $opts{infiles});
 if (scalar(@ins) != 2) { die "The number of vcf files must be two, we support merging of two files ONLY"; }
 $formats = &read_fields(\@ins);
}

### load the vcf information
my %vcf=%{load_vcf($input)};
my $imputeGT;
print STDERR scalar(keys %vcf)." entries loaded from input vcf files\n" if DEBUG;
print STDERR scalar(keys %{$vcf{calls}})." Calls loaded\n" if DEBUG;

#### these are header lines that need to be removed

our %headers2clean = (source => 1, inputs => 1, source_version => 1, content => 1, startTime => 1, cmdline => 1);

map {delete($vcf{header}{$_})} (keys %headers2clean);
print STDERR "Finished cleaning\n" if DEBUG;


my @keys = (keys %{$vcf{header}});

foreach my $head(qw/fileformat FORMAT INFO FILTER/) {
    my @newkeys = map {$_ if $_ ne $head} @keys;
    if ($vcf{header}{$head}) {
     map{$vcf{headerblock}.="##".$head."=".$vcf{header}{$head}{$_}."\n"} (sort keys %{$vcf{header}{$head}});
    }
    @keys = @newkeys;
}

FIELD:
foreach my $key(sort @keys) {
    map{if ($key eq $_){next FIELD;}} (qw/fileformat FORMAT INFO FILTER/);
    map{$vcf{headerblock}.="##".$key."=".$vcf{header}{$key}{$_}."\n"} (sort keys %{$vcf{header}{$key}});
}

my $fo = new IO::File(">$output") or die "Couldn't open file [$output] for writing";
print $fo $vcf{headerblock};
print $fo $vcf{colheaders}."\n";
map{print $fo "$vcf{calls}{$_}{line}\n"} (sort keys %{$vcf{calls}});

$fo->close;

=head2 
 this subroutine is for checking duplicate FORMAT headers
=cut

sub read_fields {

 my $vcfs = shift @_;
 my $keys = {};

 for (my $f = 0; $f < 2; $f++) {
   my $instream = $vcfs->[$f] =~/.gz$/ ? "zcat $vcfs->[$f] | " : "<$vcfs->[$f]";
   (open VCF, $instream) || die "could not open vcf file $vcfs->[$f]";
   while (<VCF>) {
   if (/^##/) {
     my($metakey,$metaval)=/^##(.*?)=(.*)/;   ### key value pairs, separated by =
     if($metaval=~/^\<.*\>$/) {
        my ($id) = $metaval =~/ID=([^,<>]+)/;
        if ($metakey eq "FORMAT") {
            $keys->{$id}++;
        }
     }
     next;
   } 
   last;
  }
  close(VCF);
 }
 return $keys;
}

=head2 load_vcf
 this subroutine will load the vcf files into a single hash
 requires both the strelka and mutect files
 strelka file will be loaded first...all records
 mutect file will be loaded second...intesecting records
 non-intersecting strelka records then removed
=cut

sub load_vcf{
	my %vcf;
	$input = shift @_;
	
                my $instream = $input =~/.gz$/ ? "zcat $input | " : "<$input";
		(open VCF, $instream) || die "could not open vcf file $input";
                my $imputeGT  = 0;
		while(<VCF>){
			chomp;
			if(/^##/){	### the headers #########################
				my($metakey,$metaval)=/^##(.*?)=(.*)/;   ### key value pairs, separated by =
                                
                                if($metaval=~/^\<.*\>$/) {
                                  my ($id) = $metaval=~/ID=([^,<>]+)/;
                                  if ($metaval=~/Description/ && $metaval!~/ID=PASS/ && $metaval!~/ID=REJECT/) {
                                        $metaval=~s/(Description=\".*?)\"/$1;source=$opts{source}\"/;
                                  } else {    ### add description key
                                      if (!/contig/ && $metaval!~/ID=PASS/ && $metaval!~/ID=REJECT/) {
                                        $metaval=~s/\>/,Description=\"source=$opts{source}\">/;
                                      }
                                  }

				  ### disabiguation code
                                  if ($metakey eq "FORMAT" && $opts{infiles}) {
                                    if ($formats->{$id} > 1) {
                                     $metaval=~s/ID=([^,<>]+)/ID=$1$opts{index}/;
                                    }
                                  } 
                                  $vcf{header}{$metakey}{$id} = $metaval;
                                } else {
                                  if ($metaval=~/Description/) {
                                        $metaval=~s/(Description=\".*?)\"/$1;source=$opts{source}\"/;
                                  } else {    ### add description key
                                        if ($metaval=~/\>/) {
                                            $metaval=~s/\>/,Description=\"source=$opts{source}\">/;
                                        } else {
                                            if ($metakey=~/germlineSnvTheta/ || $metakey=~/priorSomaticSnvRate/) {
                                                $metaval.=",Description=\"source=$opts{source}\"";
                                            }
                                        }
                                  }
                                  print STDERR "AFTER: $metaval\n" if DEBUG;
                                  $vcf{header}{$metakey}{noid}=$metaval;
                                }

			} elsif (/^#/){	### the column headerline ########################
                                my $colheader = $_;
                                $colheader =~s/FORMAT.*/FORMAT/;
                                $colheader = join("\t",($colheader, "NORMAL","TUMOR"));
				$vcf{colheaders} = $colheader;
                                if (!$vcf{header}{FORMAT}{GT}) {
                                  $vcf{header}{FORMAT}{GT} = "<ID=GT,Number=1,Type=String,Description=\"Genotype;source=$opts{source}\">";
                                  $imputeGT = 1;
                                }

			} else {### the data records ##########################
				my @f=split /\t/;
                                if ($opts{canonical} && ($f[0]=~/^chr/ && $f[0]=~/_/)) { next; } # Not clear how to deal with non-canonical chroms in model organisms (mouse?)
				my $call  = join(":",@f[0,1,3,4]);

                                ### Disambiguate FORMAT fields if we have infiles
                                if ($opts{infiles}) {
                                    my @tfs = split(":",$f[8]);
                                    for (my $tf = 0; $tf < @tfs; $tf++) {
                                       map{if($tfs[$tf] eq $_ && $formats->{$_} > 1){$tfs[$tf].=$opts{index}}} (keys %{$formats});
                                    }
                                    $f[8] = join(":",@tfs);
                                }

                                ### Calculate GT using SGT
                                if ($imputeGT) { 
                                    $f[8] = join(":",("GT",$f[8]));
                                    my($ref, $alt) = ($f[3], $f[4]);
                                    my $class=(length($ref)==1 && length($alt)==1) ? "snv" : "indel";
                                    my ($norm,$tumor) = $f[7]=~/SGT=(.*?)->(.*?);/;
                                    my ($gt_norm,$gt_tumor);
                                    
                                    if ($norm && $tumor) {
                                     if ($class eq "indel") {
                                      $gt_norm  = $norm eq "ref" ? "0/0" : $norm eq "het" ? "0/1" : "1/1";
                                      $gt_tumor = $tumor eq "ref" ? "0/0" : $tumor eq "het" ? "0/1" : "1/1";
                                     } else {
                                      $gt_norm  = "0";
                                      $gt_tumor = "."; # snv produced by MuTect are 0/1 calls 
                                     }
                                     $f[9]  = join(":", ($gt_norm, $f[9]));
                                     $f[10] = join(":", ($gt_tumor,$f[10]));
                                    }
                                    
                                }

                                $vcf{calls}{$call}{line}   = join("\t",@f);
			}	
		}
        close VCF;
	#### return the vcf record hash
	return \%vcf;
}


#! /usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use constant DEBUG=>1;

=head1 vcf_vetting.pl 
 
 This script is based on earlier work on merging mutect/strelka files
 and essentially performs the following tasks:

 * removes non-canonical contigs

 * adds GT field (dot or calculated based on SGT, if available)

 * removes tool-specific header lines

=cut

### this script will take the strelka or mutect vcf files and clean the headers, also removing non-canonical contigs
### will impute GT field if it is absent

my %opts = (canonical=>1);  ### only show canonical

### Optional filter (may be PASS to include only PASS calls)
GetOptions('filter=s'     => \$opts{filter});
my $input = shift @ARGV;

### Check if we have an input

die "Need some input to work with" if ! -e $input;

### Define canonical chromosomes
my %chroms = map{"chr".$_ => 1} (1..22,"X","Y","M");


### load the vcf information
my %vcf=%{load_vcf($input)};
my $imputeGT;
print STDERR scalar(keys %vcf)." entries loaded from input vcf files\n" if DEBUG;

#### these are headers that need to be injected
####  this should have been handled when the strelka and mutect output was filtered and annotated
our %headers2inject=(
	INFO=>{
		ANNOVAR        =>"<ID=ANNOVAR,Number=.,Type=String,Description=\"Annotation from Annovar\">",
		ANNOVAR_EXONIC =>"<ID=ANNOVAR_EXONIC,Number=.,Type=String,Description=\"Annotation from Annovar\">",
		DBSNP_GMAF     =>"<ID=DBSNP_GMAF,Number=1,Type=String,Description=\"dbSNP, global minor allele frequency\">",
		DBSNP_ALLELE   =>"<ID=DBSNP_ALLELE,Number=1,Type=String,Description=\"Allele Identified in dbSNP137\">",
		DBSNP_STRAND   =>"<ID=DBSNP_STRAND,Number=1,Type=String,Description=\"dbSNP strand\">",
		DBSNP_G5A      =>"<ID=DBSNP_G5A,Number=1,Type=String,Description=\"Identified in dbSNP G5A subset\">",
		DBSNP_G5       =>"<ID=DBSNP_G5,Number=1,Type=String,Description=\"Identified in dbSNP G5 subset\">",
		TRACK          =>"<ID=TRACK,Number=1,Type=String,Description=\"Additional tracks related to variant position\">"
	}
);

#### these are header lines that need to be removed

our %headers2clean = (source => 1, inputs => 1, source_version => 1, content => 1, startTime => 1, cmdline => 1);

foreach my $field(keys %headers2inject) {
 map {$vcf{header}{$field}->{$_} = $headers2inject{$field}->{$_}} (keys %{$headers2inject{$field}});
}
map {delete($vcf{header}{$_})} (keys %headers2clean);
print STDERR "Finished injecting and cleaning\n" if DEBUG;


#print STDERR Dumper($vcf{header});
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


print $vcf{headerblock};
print $vcf{colheaders}."\n";
map{print "$vcf{calls}{$_}{line}\n"} (sort keys %{$vcf{calls}});

#### this will load the vcf files into a single hash
#### requires both the strelka and mutect files
#### strelka file will be loaded first...all records
#### mutect file will be loaded second...intesecting records
#### non-intersecting strelka records then removed
sub load_vcf{
	my %vcf;
	$input = shift @_;
	
                my $instream = $input =~/.gz$/ ? "zcat $input | " : "<$input";
		(open VCF, $instream) || die "could not open vcf file $input";
                my $imputeGT  = 0;
		while(<VCF>){
			chomp;
			if(/^##/){   	### the headers #########################
				my($metakey,$metaval)=/^##(.*?)=(.*)/;   ### key value pairs, separated by =
                                
                                if($metaval=~/^\<.*\>$/) {
                                  my ($id)=$metaval=~/ID=([^,<>]+)/;
                                  $vcf{header}{$metakey}{$id}=$metaval;
                                } else {
                                  $vcf{header}{$metakey}{noid}=$metaval;
                                }

			} elsif (/^#/){  	### the column headerline ########################
                                my $colheader = $_;
                                $colheader =~s/FORMAT.*/FORMAT/;
                                $colheader = join("\t",($colheader, "NORMAL","TUMOR"));
				$vcf{colheaders}=$colheader;
                                if (!$vcf{header}{FORMAT}{GT}) {
                                     $vcf{header}{FORMAT}{GT} = "<ID=GT,Number=1,Type=String,Description=\"Genotype, contructed from SGT INFO via external modification\">";
                                     $imputeGT = 1;
                                }

			} else {   		#### the data records ##########################
				my @f=split /\t/;
                                if ($opts{canonical} && !$chroms{$f[0]}) { next; }
				my $call  = join(":",@f[0,1,3,4]);

                                ### Calculate GT using SGT
                                if ($imputeGT) { 
                                    $f[8] = join(":",("GT",$f[8]));
                                    my($ref, $alt) = ($f[3], $f[4]);
                                    my $class=(length($ref)==1 && length($alt)==1) ? "snv" : "indel";
                                    my ($norm,$tumor) = $f[7]=~/SGT=(.*?)->(.*?);/;
                                    my ($gt_norm,$gt_tumor);
                                    
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
				#$vcf{calls}{$call}{info}   = $f[7];
				#$vcf{calls}{$call}{format} = $f[8];
				#$vcf{calls}{$call}{normal} = $f[9];  
				#$vcf{calls}{$call}{tumor}  = $f[10];   ### 
                                
				#$vcf{calls}{$call}{main}   = join("\t",@f[0..6]);
                                $vcf{calls}{$call}{line}   = join("\t",@f);
			}	
		}
        close VCF;
	#### return the vcf record hash
	return \%vcf;
}

=head3 full_header

 function that forms a headerblock

=cut


sub full_header {
        my ($vcf,$opts)=@_;
        my ($headers)=@_;
        my @keys=sort keys %$headers;

        my $fullheader;
        ### keys to start with
        for my $key(qw/FORMAT INFO FILTER/){
                for my $id(sort keys %{$$headers{$key}}) {
                        for my $caller(qw/strelka mutect/){
                                if(my $headerline=$$headers{$key}{$id}{$caller}){
                                $fullheader.="##".$key."=".$headerline."\n";
                                }
                        }
                }
                #if(my $new=$newheaders{$key}){
                #        $fullheader.=$new;
                #}

        }

        for my $key(@keys){
                for my $id(sort keys %{$$headers{$key}}){
                        for my $caller(qw/strelka mutect/){
                                $fullheader.="##".$key."=".$$headers{$key}{$id}{$caller}."\n" if(defined $$headers{$key}{$id}{$caller});
                        }
                }
        }
        return $fullheader;
}

#!/usr/bin/env perl

################
# Perl modules #
################
use lib qw(/home/production/rtonda/Software/perl_libs/lib64/perl5/); 
use strict;
use warnings;
use Data::Dumper qw(Dumper);
use v5.10.1;	#Old perl functions allowed
use Getopt::Long;
use Template;
use Cwd;	#Get current work directory
use File::Path qw(remove_tree);	#Remove files
use File::Basename;

################################
# Figures and tables directory #
################################
my $pwd = cwd();
my $figures_folder ="$pwd/Figures";
my $tables_folder ="$pwd/Tables";
mkdir $figures_folder;
mkdir $tables_folder;
mkdir "$pwd/Report";

###################################
# Command line options processing #
###################################
my %options;
my ($inputfile, $user_fp_value,$specific_fp_value,$cosmic_db_annotation_file,$cnv_analysis,$project_title,$methods, $templatefile);
my $arg_num = scalar @ARGV;
GetOptions(\%options,
   "f=s" => \$inputfile,
   "template=s" => \$templatefile,
   "p=s", => \$user_fp_value,
   "s=s", => \$specific_fp_value,
   "cosmic=s", => \$cosmic_db_annotation_file,
   "cnv=s", =>\$cnv_analysis,
   "project=s", =>\$project_title,
   "methods:s", =>\$methods,
   "h:s",
   );
if ($options{'h'} or $arg_num < 1){usage();}
if (!$inputfile){usage();}

# Message about this program and how to use it
   sub usage
    {
        print "Script's usage description...

		usage: $0 -f file -template file [-p value] [-s value] [-project \"string\"] [-cnv] [-methods] [-cosmic file] [-h]

            -h
            	this (help) message

            -f file
            	variant call format file (.vcf) to be analyzed

            -template file
            	template Toolkit file (.tt) to be used as a template, if not defined it will use the default

            -p value
            	add extra p-values to the default p-values (1,0.05 and 0.001) that will be used for the somatic variants filtering

            -s value
            	somatic variants will be only filtered for the specified p-values on this option

            -cosmic file
            	COSMIC database file for SNPSift annotation (default \"CosmicCodingMuts_v68\")

            -cnv \"path\"
            	if the copy-number analysis was performed this options will look for its output in the specified path. If it is found it will be added to the report

            -project \"string\"
            	add the name of the project to the report title page

            -methods
            	print the methods secction as an appendix in the report (if not defined it will be not printed) 
                                
            example: $0 -f file -template file
            example: $0 -f file -template file -p value
            example: $0 -f file -template file -p value,value
            example: $0 -f file -template file -project \"string\" -cnv -methods
\n"; exit(1);
    }

#########################
# Check script's option #
#########################
#Report template file
if (!defined $templatefile){$templatefile='/project/production/DAT/apps/SOMATIC-REPORT/reporttemplate.tt';}
#Project title
if (!defined $project_title){$project_title = undef;}
#Methods appendix
if (!defined $methods){$methods = undef;}else{$methods=1;}
#CNV
my ($control_freeC_output,$cnv_output);
if (defined $cnv_analysis){
	#Search for the control-freeC output
	my @control_freeC_filename = glob $cnv_analysis . '/*.png'; #Control-freeC output path where looking for a .png file
	if (!@control_freeC_filename){die("Control-freeC output file (.png) not found.\n");}
	foreach my $out (@control_freeC_filename){
		if (-e $out) {
			$cnv_output = "$out";
		}
	}
}
else{$cnv_analysis=undef;}

#######################
# Open vcf file input #
#######################
unless (-e $inputfile) {die "No file named \"$inputfile\" found.\n";} 
if ($inputfile =~ /\.gz$/ or $inputfile =~ /\.vcf$/) {
	if ($inputfile =~ /\.gz$/){
		open(INPUT, "gunzip -c $inputfile |") or die "can't open pipe to $inputfile";
	}
	else{
		open(INPUT, "$inputfile") or die "can't open pipe to $inputfile";
	}
}

#############################
# SNPSift COSMIC Annotation #
#############################
print "SNPSift COSMIC annotation:\t";
#SNPSift has annotated chromosome M as "ChrMT", so in the case that our vcf has "ChrM" we must change it to "ChrMT"
my $filename = basename $inputfile;
my $pre_snpsift = $filename;
$pre_snpsift=~ s/.vcf.*//;
$pre_snpsift = $pre_snpsift . '.pre_snpsift' . '.vcf';

open (PRESNPSIFT, '>', "$pwd/$pre_snpsift") or die;
while (my $linedata = <INPUT>){
 	if ($linedata =~ /chrM/){;
 		my @linedata = split(/\t/, $linedata);
 		$linedata[0] =~ s/chrM/chrMT/g;
 		print PRESNPSIFT join ("\t",@linedata);
 	}
 	else{
 		print PRESNPSIFT "$linedata";
 	}
}
close($pre_snpsift);

#Now we are able to run SNPSift
my $cosmic_db;
if (defined $cosmic_db_annotation_file){$cosmic_db = $cosmic_db_annotation_file;}else{$cosmic_db="/project/production/DAT/apps/SOMATIC-REPORT/CosmicCodingMuts_v68";}
my $snpsift_output = $filename;
$snpsift_output=~ s/.vcf.*//;
$snpsift_output = $snpsift_output . '.cosm' . '.vcf';
my $runSNPSift = "java -jar /project/production/DAT/apps/SNPSIFT/1.3.4b/SnpSift.jar annotate $cosmic_db $pre_snpsift > $snpsift_output";
system ($runSNPSift) == 0 || die "system $runSNPSift failed: $?";

my $gzip= "gzip -f $snpsift_output";
system ($gzip) == 0 || die "system $runSNPSift failed: $?";

print "done\n";
#########################
# Variables declaration #
#########################
my @fp_values_array = ("1","0.05","0.0001");
my @sample_type= ("All","Control","Tumor","Shared");
my @mutation_type = ("All","SNV","INDEL");
my @effect_impact = ("High","Moderate","Low","Modifier");
my @main_statistics_array = ("AC","DP","GQ");
my @nucleotids = ("A","T","C","G");
my @outputs = ("Summary.table","Main.statistics","Zygosity","GMAF","Mutation.spectrum");

if (defined $user_fp_value){
	my @user_fp_values = split(',', $user_fp_value);
	foreach my $element(@user_fp_values){
		if ( "$user_fp_value" ~~ @fp_values_array ){
			next;
		}
		else{
			push @fp_values_array, $element;
		}
	}
}
my @all_fp_values = @fp_values_array;
my @array;
if (defined $specific_fp_value){
	@fp_values_array = split(',', $specific_fp_value);
	@array = split(',', $specific_fp_value);
	foreach my $element (@array){
		push @all_fp_values, $element;
	}
}

print "Current test p-values:",join ("\t",@fp_values_array),"\n";
print "VCF processing:\t";
#Filehandle hash
my %filehandle;
foreach my $output (@outputs){
	foreach my $fp (@fp_values_array){
		$filehandle{$output}{$fp}=undef;
	}
}

#Chromosomes
my %chrs;
foreach my $fp ( @fp_values_array ){
	$chrs{$fp}=undef;
}
#All variants summary table hash
my %mutations;
foreach my $element (@sample_type ) {
	foreach (@mutation_type){
		$mutations{$element}{$_}=0;
	}
}

my %impact;
foreach my $element ( @sample_type ) {
	foreach (@effect_impact){
		$impact{$element}{$_}=0;
	}
}

#Somatic variants summary table hash
my %fp_hash; 
foreach my $fp ( @fp_values_array ) {
	foreach my $sample (@sample_type){
		foreach (@mutation_type){
			$fp_hash{$fp}{$sample}{$_}=0;
		}
	}
}

my %fp_hash_impact; 
foreach my $fp ( @fp_values_array ) {
	foreach my $sample (@sample_type){
		foreach (@effect_impact){
			$fp_hash_impact{$fp}{$sample}{$_}=0;
		}
	}
}

#Somatic variants candidates hash
my %cosmic;
%cosmic= map { $_ => 0 } @all_fp_values;
my %somatic_candidates = %cosmic;

#Main statistic hash
my %main_statistics_hash; 
foreach my $fp ( @fp_values_array ) {
	foreach my $statistic (@main_statistics_array){
		$main_statistics_hash{$fp}{$statistic}{Control} = [undef];
		$main_statistics_hash{$fp}{$statistic}{Tumor} = [undef];
	}
}

#Mutation spectrum hash
my %mutation_spectrum_hash;
foreach my $fp (@fp_values_array){
	foreach my $r (@nucleotids){
		foreach my $a (@nucleotids){
			if ($r ne $a){
				my $SNP= "$r".":"."$a";
				$mutation_spectrum_hash{$fp}{$SNP}=0;
			}
		}
	}
	$mutation_spectrum_hash{$fp}{INDELs}=0
}

my %zygosity;
my @sample=("Control","Tumor");
my @genotype_array=("0/0","0/1","1/1");
foreach my $fp (@fp_values_array){
	foreach my $sample (@sample){
		foreach my $genotype (@genotype_array){
			$zygosity{$fp}{$sample}{$genotype}=0;
		}
	}
}

my %outputs;
###########
# OUTPUTS #
###########
foreach my $fp (@fp_values_array){
	my $fp_summary_table = "summarytable"."$fp".".txt";
	my $main_statistics= "mainstatistics"."$fp".".txt";
	my $mutation_spectrum = "mutation.spectrum."."$fp".".txt";
	my $gmaf_table = "gmaf."."$fp".".txt";
	my $zygosity_table = "zygosity."."$fp".".txt";
	my $variants_filtered = "$snpsift_output".'.fp'."$fp".'.vcf';
	$outputs{$fp} = [ $fp_summary_table, $main_statistics, $mutation_spectrum, $gmaf_table, $zygosity_table, $variants_filtered];
	open ($filehandle{Summary_table}{$fp}, ">","$tables_folder/$fp_summary_table") or die "Can not open $fp_summary_table";
	open ($filehandle{Main_statistics}{$fp}, ">", "$tables_folder/$main_statistics") or die "Can not open $main_statistics";
	open ($filehandle{Zygosity}{$fp}, ">", "$tables_folder/$zygosity_table") or die "Can not open $zygosity_table";
	open ($filehandle{GMAF}{$fp}, ">", "$tables_folder/$gmaf_table") or die "Can not open $gmaf_table";
	open ($filehandle{Mutation_spectrum}{$fp}, ">", "$tables_folder/$mutation_spectrum") or die "Can not open $main_statistics";
	open ($filehandle{Variants_filtered}{$fp}, ">", "$tables_folder/$variants_filtered") or die "Can not open $variants_filtered";
}
#Specific header from each output table
foreach my $element (@fp_values_array){
	print {$filehandle{Main_statistics}{$element}}"Chr\tPosition\tControl GT\tTumor GT\tControl VAF\tTumor VAF\tControl DP\tTumor DP\tFP value\tControl GQ\tTumor GQ\tControl AC\tTumor AC\n";
	print {$filehandle{Zygosity}{$element}} "Genotype\tControl\tTumor\n";
	print {$filehandle{GMAF}{$element}} "Chr\tGMAF\n";
}

##################
# vcf processing #
##################
open (INPUT,"gunzip -c $snpsift_output |") or die "can't open pipe to $snpsift_output";
while (my $linedata = <INPUT>) {
	chomp $linedata;

 	if ($linedata =~ /chrMT/){;$linedata =~ s/chrMT/chrM/g;}

 	#Print meta-info and header lines in the vcf outputs from the variants filtered by each p-value.
	if ($linedata =~ /^#/){
		foreach my $fp_value (@fp_values_array){
			print {$filehandle{Variants_filtered}{$fp_value}}"$linedata\n";
		}
	}

	else{

		#Vcf fields processing
		my ($chr, $pos, $id, $refal, $altal, $quality, $filter, $info, $format, $controlformat, $tumorformat) = split(/\t/, $linedata);
		my @infoarray = split (';', $info);
		my @formatarray = split(':', $format);
		my @idarray = split(';', $id);
		
		#Filters
		if ($filter !~ /mrd10/){	#reached a minimum read depth per sample of 10 reads
			#Genotype				
			my ($indexGT) = grep { $formatarray[$_] =~ /GT/ } 0..$#formatarray;	
			my @control = split(':', $controlformat);
			my @tumor= split(':', $tumorformat);
			my $GTcontrol = $control[$indexGT];
			my $GTtumor = $tumor[$indexGT];

			# #All mutations count
			MutationTypeCounter(\%mutations,\$info,"All");
			if ($GTcontrol !~ /0\/0/) {MutationTypeCounter(\%mutations,\$info,"Control");}
			if ($GTtumor !~ /0\/0/) {MutationTypeCounter(\%mutations,\$info,"Tumor");}
			if ($GTcontrol !~ /0\/0/ and $GTtumor !~ /0\/0/) {MutationTypeCounter(\%mutations,\$info,"Shared");}

			#Effect impact
			ImpactCounter(\%impact,\$info,"All");
			if ($GTcontrol !~ /0\/0/) {ImpactCounter(\%impact,\$info,"Control");}
			if ($GTtumor !~ /0\/0/) {ImpactCounter(\%impact,\$info,"Tumor");}
			if ($GTcontrol !~ /0\/0/ and $GTtumor !~ /0\/0/) {ImpactCounter(\%impact,\$info,"Shared");}
			
			#Somatic filter for somatic variants only
			if ($info =~ /FP=/){
				my ($fp_value) = grep(/FP=/,@infoarray);
				$fp_value=~ s/FP=//;

				#Somatic candidates
				foreach my $key (keys %somatic_candidates) {
					if ($fp_value <= $key){
						$somatic_candidates{$key}++;
						#Cosmic annotation
						if ($id =~ /COSM/){		
						$cosmic{$key}++;
						}
					}
				}

				#Analysis for each fp.value
				foreach my $fp_test_value ( @fp_values_array ) {
					#FILTER FOR SOMATIC VARIANTS (FP and mrd10)
					if($fp_value <= $fp_test_value){ 
						if (!defined $chrs{$fp_test_value}{$chr}){$chrs{$fp_test_value}{$chr}=0;}else {$chrs{$fp_test_value}{$chr}++;}	#Push chromosomes names to a hash
 					
 					#Print variants data in the vcf outputs from the variants filtered by each p-value.
 					print {$filehandle{Variants_filtered}{$fp_test_value}}"$linedata\n";

						#All mutations count
						MutationTypeCounter_fpvalue(\%fp_hash,$fp_test_value,\$info,"All");
						if ($GTcontrol !~ /0\/0/) {MutationTypeCounter_fpvalue(\%fp_hash,$fp_test_value,\$info,"Control");}
						if ($GTtumor !~ /0\/0/) {MutationTypeCounter_fpvalue(\%fp_hash,$fp_test_value,\$info,"Tumor");}
						if ($GTcontrol !~ /0\/0/ and $GTtumor !~ /0\/0/) {MutationTypeCounter_fpvalue(\%fp_hash,$fp_test_value,\$info,"Shared");}

						#Effect impact
						ImpactCounter_fpvalue(\%fp_hash_impact,$fp_test_value,\$info,"All");
						if ($GTcontrol !~ /0\/0/) {ImpactCounter_fpvalue(\%fp_hash_impact,$fp_test_value,\$info,"Control");}
						if ($GTtumor !~ /0\/0/) {ImpactCounter_fpvalue(\%fp_hash_impact,$fp_test_value,\$info,"Tumor");}
						if ($GTcontrol !~ /0\/0/ and $GTtumor !~ /0\/0/) {ImpactCounter_fpvalue(\%fp_hash_impact,$fp_test_value,\$info,"Shared");}

						#Main statistics
						foreach my $statistic (@main_statistics_array){
							my ($index) = grep { $formatarray[$_] =~ /$statistic/ } 0..$#formatarray;
							my $statistic_control = $control[$index];
							my $statistic_tumor = $tumor[$index];
							@{ $main_statistics_hash{$fp_test_value}{$statistic}{Control} } = $statistic_control;
							@{ $main_statistics_hash{$fp_test_value}{$statistic}{Tumor} } = $statistic_tumor;
							# push @{ $main_statistics_hash{$fp_test_value}{$statistic}{Tumor} }, $statistic_tumor;
						
							my ($VAF_control, $DP_control, $VAF_tumor, $DP_tumor);
							if ($statistic =~ /AC/){
								my ($reference_control, $alternative_control) = split ",",$statistic_control;
							    if ($reference_control == 0 and $alternative_control == 0){
							    	$DP_control=0;
							    	$VAF_control= 0;
							    }
							    else{
							    	$DP_control = $reference_control+$alternative_control;
							    	$VAF_control = $alternative_control/($reference_control+$alternative_control);
							    }
							    $VAF_control = eval sprintf('%.3f', $VAF_control);

							    my ($reference_tumor, $alternative_tumor) = split ",",$statistic_tumor;
							    if ($reference_tumor == 0 and $alternative_tumor == 0){
							    	$DP_tumor=0;
							    	$VAF_tumor= 0;
							    }
							    else{
							    	$DP_tumor = $reference_tumor+$alternative_tumor;
							    	$VAF_tumor = $alternative_tumor/($reference_tumor+$alternative_tumor);
							    }
							    $VAF_tumor = eval sprintf('%.3f', $VAF_tumor);

							    print {$filehandle{Main_statistics}{$fp_test_value}}"$chr\t$pos\t$GTcontrol\t$GTtumor\t$VAF_control\t$VAF_tumor\t$DP_control\t$DP_tumor\t$fp_value";
							}

							if ($statistic =~ /GQ/){print {$filehandle{Main_statistics}{$fp_test_value}}"\t$statistic_control\t$statistic_tumor\t@{ $main_statistics_hash{$fp_test_value}{AC}{Control} }\t@{ $main_statistics_hash{$fp_test_value}{AC}{Tumor}}\n";}

						}
						#GMAF
						if ($id =~ /GMAF=/ ){
							my ($gmaf_value) = grep(/GMAF=/,@idarray);
							$gmaf_value=~ s/GMAF=//;
							print {$filehandle{GMAF}{$fp_test_value}} "$chr\t$gmaf_value\n";
					 	}

					 	#Mutation spectrum
					 	MutationSpectrum(\%mutation_spectrum_hash, $fp_test_value,$refal,$altal, $GTcontrol, $GTtumor);
					 	if ($info =~ /INDEL/){$mutation_spectrum_hash{$fp_test_value}{INDELs}++}

					 	#Zygosity
					 	for ($GTcontrol) {
						when (/0\/0/) { $zygosity{$fp_test_value}{Control}{'0/0'}++ }
						when (/0\/1/) { $zygosity{$fp_test_value}{Control}{'0/1'}++  }
						when (/1\/1/) { $zygosity{$fp_test_value}{Control}{'1/1'}++  }
						}
						for ($GTtumor) {
						when (/0\/0/) { $zygosity{$fp_test_value}{Tumor}{'0/0'}++  }
						when (/0\/1/) { $zygosity{$fp_test_value}{Tumor}{'0/1'}++ }
						when (/1\/1/) { $zygosity{$fp_test_value}{Tumor}{'1/1'}++ }
						}
					}
				}
			}
		}
	}
}

###############################
# All mutations summary table #
###############################
my $all_summary_table = "$tables_folder/all.summary.table.txt";
open (SUMMARYTABLE,">", "$all_summary_table") or die ("Couldn't open all.summary.table.txt");
#Mutation type
print SUMMARYTABLE "Type\tAll\tControl\tTumor\tShared\n";
foreach my $mutation (@mutation_type){
	print SUMMARYTABLE "$mutation";
	print SUMMARYTABLE "\t$mutations{All}{$mutation}\t$mutations{Control}{$mutation}\t$mutations{Tumor}{$mutation}\t$mutations{Shared}{$mutation}\n";
}
#Effect impact
print SUMMARYTABLE "Effectimpact\t-\t-\t-\t-\n";
foreach my $effect (@effect_impact){
	print SUMMARYTABLE "$effect\t$impact{All}{$effect}\t$impact{Control}{$effect}\t$impact{Tumor}{$effect}\t$impact{Shared}{$effect}\n";
}
print SUMMARYTABLE"\n";

###########################################
# Somatic variants candidates by FP value #
###########################################
my $somatic_candidates_table = "$tables_folder/FPthreshold.table.txt";
open (SOMATICTHRESHOLD,">", "$somatic_candidates_table") or die ("Couldn't open FPthreshold.table.txt");
print SOMATICTHRESHOLD "p-value\tVariants\tdbCOSMIC\n";
for (keys %somatic_candidates){
	print SOMATICTHRESHOLD "$_\t$somatic_candidates{$_}\t$cosmic{$_}\n"
}

#####################################################
# Somatic mutations summary table for each fp-value #
#####################################################
#Mutation type
foreach my $element (@fp_values_array){
	print {$filehandle{Summary_table}{$element} }"Type\tAll\tControl\tTumor\tShared\n";
}

foreach my $key (keys %fp_hash){
	foreach my $mutation (@mutation_type){
		print {$filehandle{Summary_table}{$key}} "$mutation";
		foreach my $sample (@sample_type){
			print {$filehandle{Summary_table}{$key}} "\t$fp_hash{$key}{$sample}{$mutation}";
		}
		print {$filehandle{Summary_table}{$key}} "\n";
	}
}
#Effect impact
foreach my $test_value (@fp_values_array){
	print {$filehandle{Summary_table}{$test_value}} "Effectimpact\t-\t-\t-\t-\n";
} 
foreach my $key (keys %fp_hash_impact){
	foreach my $effect (@effect_impact){
		print {$filehandle{Summary_table}{$key}} "$effect";
		foreach my $sample (@sample_type){
			print {$filehandle{Summary_table}{$key}} "\t$fp_hash_impact{$key}{$sample}{$effect}";
		}
		print {$filehandle{Summary_table}{$key}} "\n";
	}

}

#####################
# Mutation spectrum #
#####################
foreach my $fp (@fp_values_array){
	print {$filehandle{Mutation_spectrum}{$fp}} "Type\tNumber\n";
	foreach my $key (keys %{$mutation_spectrum_hash{$fp}}){
		print {$filehandle{Mutation_spectrum}{$fp}} "$key\t$mutation_spectrum_hash{$fp}{$key}\n";
	}
}

########
# GMAF #
########
# foreach my $fp_value_array (@fp_values_array){
# 	my $count=0;
# 	close ($filehandle{GMAF}{$fp_value_array});
# 	my $file1 = "$tables_folder/gmaf."."$fp_value_array".".txt";
# 	open(my $file, "<", $file1);
# 	while (my $linedata = <$file>){
# 		$count++;
# 	}
# 	if ($count < 2){
# 		open (my $gmaf_output, ">", "$file1");
# 		foreach my $keys (keys $chrs{$fp_value_array}){
# 			print $gmaf_output "Chr\tGMAF\n$keys\t0\n";
# 		}
# 	}
# }

############
# Zygosity #
############
my @zygosity_array = ("ref/ref","ref/alt","alt/alt");

foreach my $fp (@fp_values_array){
	foreach my $element (@zygosity_array){
		print {$filehandle{Zygosity}{$fp}} "$element";
		if ($element eq "ref/ref"){
			print{$filehandle{Zygosity}{$fp}} "\t$zygosity{$fp}{Control}{'0/0'}\t$zygosity{$fp}{Tumor}{'0/0'}\n";
		}
		if ($element eq "ref/alt"){
			print {$filehandle{Zygosity}{$fp}}"\t$zygosity{$fp}{Control}{'0/1'}\t$zygosity{$fp}{Tumor}{'0/1'}\n";
		}
		if ($element eq "alt/alt"){
			print {$filehandle{Zygosity}{$fp}}"\t$zygosity{$fp}{Control}{'1/1'}\t$zygosity{$fp}{Tumor}{'1/1'}\n";
		}
	}
	print {$filehandle{Zygosity}{$fp}}"\n";
}

print "done\n";
###############
# Subroutines #
###############
sub MutationTypeCounter {
	my ($hashref,$inforef,$tag)=@_;

	for ($$inforef) {
		$$hashref{$tag}{"All"}++;
		when ($$inforef =~ /INDEL/) {$$hashref{$tag}{"INDEL"}++}
		when ($$inforef !~ /INDEL/) {$$hashref{$tag}{"SNV"}++}
	}
	return \$hashref;
}

sub ImpactCounter {
	my ($hashref,$inforef,$tag)=@_;

	for ($$inforef) {
		when ($$inforef =~ /HIGH/) {$$hashref{$tag}{"High"}++}
		when ($$inforef =~ /MODERATE/) {$$hashref{$tag}{"Moderate"}++}
		when ($$inforef =~ /LOW/) {$$hashref{$tag}{"Low"}++}
		when ($$inforef =~ /MODIFIER/) {$$hashref{$tag}{"Modifier"}++}
	}
	return \$hashref;
}

sub MutationTypeCounter_fpvalue {
	my ($fp_hash,$fpvalue,$inforef,$tag)=@_;

	for ($$inforef) {
		$$fp_hash{$fpvalue}{$tag}{"All"}++;
		when ($$inforef =~ /INDEL/) {$$fp_hash{$fpvalue}{$tag}{"INDEL"}++}
		when ($$inforef !~ /INDEL/) {$$fp_hash{$fpvalue}{$tag}{"SNV"}++}
	}
	return $fp_hash;
}

sub ImpactCounter_fpvalue {
	my ($fp_hash_impact,$fpvalue,$inforef,$tag)=@_;

	for ($$inforef) {
		when ($$inforef =~ /HIGH/) {$$fp_hash_impact{$fpvalue}{$tag}{"High"}++}
		when ($$inforef =~ /MODERATE/) {$$fp_hash_impact{$fpvalue}{$tag}{"Moderate"}++}
		when ($$inforef =~ /LOW/) {$$fp_hash_impact{$fpvalue}{$tag}{"Low"}++}
		when ($$inforef =~ /MODIFIER/) {$$fp_hash_impact{$fpvalue}{$tag}{"Modifier"}++}
	}
	return $fp_hash_impact;
}

sub MutationSpectrum {
 	my ($hashref, $fpvalue, $r, $a, $GTcontrol, $GTtumor)=@_;
 	if (length($r) == 1 && length($a) == 1){
 		my $SNP = "$r".":"."$a";
 		if($GTcontrol eq "0/1" and $GTtumor eq "0/0"){$SNP = "$a".":"."$r";}
 		if ($GTcontrol eq "1/1" and $GTtumor eq "0/0"){$SNP = "$a".":"."$r";}
 		if ($GTcontrol eq "1/1" and $GTtumor eq "0/1"){$SNP = "$a".":"."$r";}
	 	if(defined $$hashref{$fpvalue}{$SNP}){
	 		{$$hashref{$fpvalue}{$SNP}++};
	 	}
 	}
 	if (length($r) == 1 && $a =~ /[^,]+(,[^,]+)+/g){	#extract only the comma-separated strings
 		my @alt_array = split (',',$a);
 		foreach my $element (@alt_array){
 			my $SNP = "$r".":"."$element";
 			if($GTcontrol eq "0/1" and $GTtumor eq "0/0"){$SNP = "$a".":"."$r";}
 			if ($GTcontrol eq "1/1" and $GTtumor eq "0/0"){$SNP = "$element".":"."$r";}
 			if ($GTcontrol eq "1/1" and $GTtumor eq "0/1"){$SNP = "$element".":"."$r";}
 			if(defined $$hashref{$fpvalue}{$SNP}){
 				{$$hashref{$fpvalue}{$SNP}++};
 			}
 		}
 	}
}

############
# Template #
############
my $config = {
#    INCLUDE_PATH => '/search/path',  # or list ref
     INTERPOLATE  => 0,               # expand "$var" in plain text
     POST_CHOMP   => 1,               # cleanup whitespace
     ABSOLUTE     => 1, 
#    PRE_PROCESS  => 'header',        # prefix each template
#    EVAL_PERL    => 1,               # evaluate Perl code blocks
};

# create Template object
my $template = Template->new($config);

#Run rainfallscript
my $rainfall_plot = "$snpsift_output";
my $rainfall_raul_path = '/project/production/DAT/apps/SOMATIC-REPORT/Rainfall/Rainfall.R';
my $hash_ref = \%outputs;

my $data={
	title=>"Somatic mutations report",
	logo=>"Logo/logo.jpeg",
	project=>$project_title,
	allsummarytable=>$all_summary_table,
	somaticcandidates=>$somatic_candidates_table,
	rainfallplot=> $rainfall_plot,
	rainfallscript=>$rainfall_raul_path,
	figuresfolder=>$figures_folder,
	fpvalues=>\@fp_values_array,
	fpvalueshash=> $hash_ref,
	cnv_analysis=>$cnv_output,
	methods=>$methods,
};

# specify input filename, or file handle, text reference, etc.
unless (-e $templatefile) {die "No file named \"$templatefile\" found, check if your template (.tt) file exists please!\n";}
my $templatefilename = basename $templatefile;


# process input template, substituting variables
my $rwn_file= $templatefilename;
$rwn_file =~ s/\.[^.]*$//;
$rwn_file = $rwn_file . ".Rnw";
#print "\n$rwn_file\n\n";
$template->process($templatefile, $data, $rwn_file)
    || die $template->error("Template error\n");

my $RSweave="Rscript -e 'Sweave(\"$rwn_file\")'";
system($RSweave);

my $texfile = $rwn_file;
$texfile =~ s/\.[^.]*$//;
$texfile = $texfile . ".tex";
my $pdflatex="pdflatex $texfile";
system($pdflatex);

$texfile =~ s/\.[^.]*$//;
my $bibtex = "bibtex $texfile";
system($bibtex);

$texfile = $texfile . ".tex";
my $pdflatex2="pdflatex $texfile";
system($pdflatex2);

my $pdflatex3="pdflatex $texfile";
system($pdflatex3);

$texfile =~ s/\.[^.]*$//;
my $pdffile = $texfile . ".pdf";
my $copy_report="cp $pdffile Report";
system($copy_report);


#Delete needless files once report is generated
my $delete_temp='rm *.aux *.out *.toc *.log *.lof *.lot *.Rnw *.tex *.bbl *.blg *.pdf *.imd.mut.df.gz';
system($delete_temp);
$snpsift_output="$snpsift_output".'.gz';
remove_tree("$snpsift_output","$pre_snpsift");

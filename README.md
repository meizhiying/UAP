# UAP
Duplex UMI split and analysis package

## About UAP  

Name: UMI Analysis Package(UAP)  
Version: 1.1.0  
Function:  

	a.Split Duplex UMI with PE fastq files.  
	b.Calls consensus sequences from reads with the same unique molecular tag and applys sequencing error correction with unique molecular tag.  
	c.Calculates error rate from aligned bam file.  

Environment:linux  
Help information:run UAP in command line will get detailed information  

## Installation  

To clone the repository:  

	git clone https://github.com/meizhiying/UAP.git

Installation command line:  

	sh path_to_UAP/custom_install.sh

Set the UAP_HOME environment variable to UAP path

	export UAP_HOME=path_to_UAP

Set the JAVA_HOME environment variable to Java 8 path

	export JAVA_HOME=path_to_java8

## Params of UAP  

Run path_to_UAP/bin/UAP will get the help information  

    UMI analysis package.
    ---------------------
    Version:1.1.0

    Usage:
            UAP <command> [arguments]

    Command list:
            AnnoFastqWithUMI <Annotates fastq file with UMI information>
                    Arguments:
                    -f <FILE>       Input fastq1 file
                    -r <FILE>       Input fastq2 file
                    -e <INT>        Maximum mismatch bases in UMIs[Default:0]
                    -o <STR>        Output prefix
                    -d <STR>        Output directory

            BamConsensusAnalysis <Annotates existing BAM files and filter consensus reads>
                    Arguments:
                    -i <FILE>       Input bam file
                    -o <STR>        Output prefix
                    -d <STR>        Output directory
                    -b <FILE>       Target bed file
                    -m <INT>        Minimum mapping quality
                    -s <INT>        Family size threshold
                    -u <STR>        UMI type [Duplex/Single]
                    -B <FLOAT>      Maximum base error rate in UMI family
                    -F <STR>        Reference path
                    -r <INT>        [0/1] 1 means remove all Temp files. Default:1

            ErrorRateStats <Calculates the error rate by read position on coordinate sorted mapped BAMs>
                    Arguments:
                    -i <FILE>       Input bam file
                    -v <VCF>        Optional VCF file of variant sites to ignore
                    -F <STR>        Reference path
                    -m <INT>        Minimum mapping quality
                    -o <STR>        Output prefix

## Requirements  

Before running UAP, you need to make sure that several pieces of software  
and/or modules are installed on the system:  

* Java 8 is required. We recommend either of the following: 
	* OpenJDK 8 with Hotspot from [AdoptOpenJdk](https://adoptopenjdk.net/)
	* [OracleJDK 8](https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
which requires an Oracle account to download and comes with restrictive [license conditions](https://www.oracle.com/downloads/licenses/javase-license1.html).
* GNU make 
* C compiler (e.g. gcc or clang) 
* R & ggplot2  


## Running the UAP  

UAP has three modules: AnnoFastqWithUMI, BamConsensusAnalysis, ErrorRateStats.The read id in BamConsensusAnalysis's input bam and AnnoFastqWithUMI's output bam must be in the same format.  

Test data are saved in path_to_UAP/test  

	Test_1.fq.gz: read1 fastq file from PE sequencing data  
	Test_2.fq.gz: read2 fastq file from PE sequencing data  
	Test.bam: Test bam file  
	Test.bed: Test bed file  


1.AnnoFastqWithUMI:  

This module will complete UMI split and annotation, command line:

	UAP AnnoFastqWithUMI -f Test_1.fq.gz -r Test_2.fq.gz -o Test -d Split_test -e 0

Three files named Test_1.splitUMI.fq.gz Test_1.splitUMI.fq.gz and Test.umi.log will be generated.  
Output fastq file contains reads with UMI information,Test.umi.log file contains UMI split rate  

2.BamConsensusAnalysis:  

This module will apply sequencing error correction with unique molecular tag, command line:  

	UAP BamConsensusAnalysis -i Test.bam -o Test -d Bam_test -b test.bed -F hg19.fa -m 10 -s 3 -u Duplex -B 0.1 -r 1

Three files named Test.target.consensus.bam Test.target.consensus.bai and Test.family_size.xls will be generated.  
Test.target.consensus.bam contains reads generated from UMI corection,Test.family_size.xls file contains UMI family size distribution.  

3.ErrorRateStats  

This module will calculate sequencing error rate from bam file, command line:  

	UAP ErrorRateStats -i Test.bam -o Error_stat_test/Test -F hg19.fa -m 10

Two files named Test.error_rate_by_read_position.txt and Test.error_rate_by_read_position.pdf will be generated  
Test.error_rate_by_read_position.txt contains detail error rate statistics,Test.error_rate_by_read_position.pdf contains plots base on Test.error_rate_by_read_position.txt  

#!/bin/bash

set -e

workdir=`pwd`
toolsdir=$workdir/tools
mkdir -p $toolsdir

BWA_URL="http://jaist.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2"
PICARD_URL="https://github.com/broadinstitute/picard/releases/download/2.21.6/picard.jar"
SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2"
FGBIO_URL="https://github.com/fulcrumgenomics/fgbio/releases/download/1.0.0/fgbio-1.0.0.jar"

DOWNLOAD_URL()
{
	local URL=$1
	local name=$2
	local unzip_dir=$3
	local unzip_info=$4

	if [ -f "$toolsdir/$name" ]
	then
		rm "$toolsdir/$name"
	fi

	echo "Downloading: $URL"

	wget -c -P $toolsdir $URL

	echo "Downloading completed successfully"

	if [[ ${unzip_info} != "F" ]]
	then
		rm -rf $toolsdir/${unzip_dir}
		tar jxf $toolsdir/$name -C $toolsdir
	fi
}

DOWNLOAD_URL ${BWA_URL} "bwa-0.7.17.tar.bz2" "bwa-0.7.17" "T"
DOWNLOAD_URL ${PICARD_URL} "picard.jar" "picard" "F"
DOWNLOAD_URL ${SAMTOOLS_URL} "samtools-1.10.tar.bz2" "samtools-1.10" "T"
DOWNLOAD_URL ${FGBIO_URL} "fgbio-1.0.0.jar" "fgbio" "F"

echo "Installing bwa ..."
cd $toolsdir/bwa-0.7.17
make
if [ -f $toolsdir/bwa ]
then
	rm $toolsdir/bwa
fi
cd .. && ln -s bwa-0.7.17/bwa
echo "Installing bwa successfully"

echo "Installing samtools ..."
cd $toolsdir/samtools-1.10
./configure --prefix=$toolsdir/samtools-1.10 --disable-lzma
make
make install
if [ -f $toolsdir/samtools ]
then
	rm $toolsdir/samtools
fi
cd .. && ln -s samtools-1.10/samtools
echo "Installing samtools successfully"

chmod 755 $toolsdir/picard.jar
chmod 755 $toolsdir/fgbio-1.0.0.jar

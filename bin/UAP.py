#!/usr/bin/env python3

import os
import sys
import re
import regex
import gzip
import time
import pysam
from subprocess import check_call
from multiprocessing import Pool
from optparse import OptionParser

class AnnoFastqWithUMI(object):
	def __init__(self):
		if not optlist.opt_f:
			print('==ERROR: missing -f arguments==')
			usage()
		if not optlist.opt_r:
			print('==ERROR: missing -r arguments==')
			usage()
		if not optlist.opt_o:
			print('==ERROR: missing -o arguments==')
			usage()
		if not optlist.opt_d:
			print('==ERROR: missing -d argements==')
			usage()
		if not os.path.exists(optlist.opt_d):
			os.makedirs(optlist.opt_d)
		self.fq1 = optlist.opt_f
		self.fq2 = optlist.opt_r
		self.outprefix = optlist.opt_o
		self.outdir = optlist.opt_d
		if optlist.opt_e:
			self.err = optlist.opt_e
		else:
			self.err = 0
	def check_umi(self,umi_list,seq,err):
		result = 'NA'
		if seq in umi_list:
			result = seq
		if result == 'NA' and err > 0:
			for umi in umi_list:
				check_list = regex.findall('(%s){e<=%s}'%(seq,err),umi)
				if len(check_list) > 0:
					result = umi
					break
		return result
	def gzip_fq(self,fastq,log):
		file_log = open(log,'a')
		file_log.write('[ '+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+' ] '+'Gzip %s (%s) ...\n'%(fastq, os.getpid()))
		start = time.time()
		check_call('gzip %s'%(fastq),shell=True)
		end = time.time()
		file_log.write('[ '+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+' ] '+'Gzip %s runs %0.2f mins.\n'%(fastq, (end - start)/60))
		file_log.close()
		return
	def annofq(self):
		umi_list = ['TCTCG', 'ACGGA', 'GCATT', 'TGGCA', 'CGGTA', 'CTAGG', 'CTTGA', 'AGGTG', 'TGGAT', 'CCATT', 'TCCAG', 'TACGT', 'AAGCT', 'AGTAT', 'ACTTC', 'CGTAG', 'AGTCG', 'ATATG', 'CTATG', 'ATGAT', 'AGAGC', 'GTGAT', 'TCAAT', 'TCCAC', 'CTCGT', 'TCCAT', 'TTGAG', 'GACAG', 'TGCCA', 'GTTCA', 'CATCC', 'AAGTT', 'CGAGA', 'AATTC', 'CGAGT', 'AATTG', 'TTGAA', 'CCGTT', 'TTAGC', 'TTACA', 'GGCAT', 'TTAAC', 'TTACC', 'AATCC', 'CCGTA', 'TTAGG', 'GGCTT', 'TTGAC', 'CCAGT', 'AATCG', 'TTCAA', 'AATGT', 'GGACA', 'AATCT', 'TTAGA', 'TTCAC', 'CCATG', 'TTGCA', 'TTCCA', 'AATGG', 'AATGC', 'GGCAA', 'GGATT', 'GTCAA', 'TGTCG', 'CAGTT', 'GTCCA', 'CATAC', 'TGCGT', 'TGGTG', 'CATCA', 'GTACG', 'GTCGT', 'GGTGT', 'GTCGA', 'TGGTC', 'CATCG', 'ACGTT', 'TGTGG', 'GAGTC', 'TCCAA', 'ACCGA', 'ATACG', 'CAACT', 'CTAAG', 'TACAC', 'CTGTC', 'GTAAT', 'AGCGT', 'ACTGG', 'GCACA', 'ATGGT', 'GACCT', 'AGGCA', 'CTCCT', 'TACCG', 'TAGTG', 'AGATT', 'TCTAC', 'CACGT', 'GGTAC', 'TAGGT', 'TCACA', 'GATGC', 'TCGAA', 'GAGTT', 'ATCTA', 'CTTAA', 'CTGAA', 'GCTGA', 'CTAGC', 'GCGAA', 'CTCTG', 'GCAGA', 'AGGAC', 'CTGAG', 'CGATG', 'GTCTC', 'CATAG', 'ACGAC', 'TGGTA', 'TGCTC', 'CATGT', 'GTACA', 'CATAT', 'ACACC', 'CATTA', 'TGAGC', 'ACGAG', 'GTCCT', 'CATGA', 'GTCTA', 'TGAGT', 'TGCGA', 'ACGAT', 'CATTC', 'CATTG', 'TGATC', 'TGAGG', 'ACGCA', 'GTAAC', 'GTCAT', 'TGCTG', 'GTAAG', 'CATGC', 'TGCTA', 'ACGCT', 'GTATA', 'ACACG', 'CATGG', 'CATCT', 'GTATC', 'GTCAC', 'GTCTG', 'CAGAT', 'GGATG', 'GTACT', 'CAGGA', 'TAATG', 'CCTTG', 'GTAGA', 'ACCAC', 'GTATG', 'ACACT', 'TGTGA', 'TGTCA', 'TGCCT', 'AACTG', 'GCTAC', 'GACAC', 'CACAG', 'CTCTA', 'ATCCT', 'TAGAG', 'TACGC', 'GAATG', 'TGTAA', 'ATGTA', 'ATCCG', 'AGTCT', 'TGTCC', 'TCTCA', 'TATAC', 'TCAAC', 'GGTGA', 'ACGTG', 'ATGTG', 'ACGTA', 'ATGCT', 'CTACG', 'CTAAT', 'ATGCG', 'GCTTG', 'ACTTG', 'GTTAA', 'ACCAG', 'CAGTC', 'ACGGT', 'CGAAC', 'GCTGT', 'CAGTG', 'ACAGG', 'CAGTA', 'TGGCT', 'CACTT', 'ACAGT', 'TGTAC', 'GTTAC', 'CAGGT', 'GTAGT', 'TGACC', 'TGATG', 'ACCAT', 'TGCAG', 'TGACA', 'ACTAG', 'TGTGC', 'GTAGC', 'CAGCT', 'GTTAG', 'GTTGA', 'ACATC', 'ACAGC', 'CGCAA', 'GATGT', 'ACATG', 'ATGGC', 'GCTTA', 'TAGCA', 'CGAAG', 'ACATT', 'CGAAT', 'GATTG', 'CTGCT', 'TCAGT', 'AGCTT', 'TAGTA', 'ATTCC', 'CGATA', 'CAATC', 'TCCTA', 'GGCTA', 'TATCG', 'TAATC', 'ACTCG', 'TACTG', 'AGACG', 'GAGAG', 'GTGCT', 'GTGCA', 'GGACT', 'GGTAT', 'GCTTC', 'GCTAG', 'AGACT', 'GTTGG', 'CGCTA', 'ATACC', 'AGCGA', 'CGGAT', 'ATCGC', 'CCTGT', 'CTCGA', 'TAGAT', 'CTCAC', 'TTGCC', 'GATTC', 'CGACA', 'GATAT', 'CGACT', 'ATGTC', 'TTCAG', 'GATAC', 'GATTA', 'GATCT', 'TAGAC', 'GCTCT', 'ACCGT', 'GCTCA', 'GATCG', 'ACCTC', 'GATGG', 'TACCA', 'ATTAC', 'CGATT', 'ACCTT', 'GATGA', 'TAGCT', 'ATTAG', 'TTCGC', 'ACCTG', 'CGCAT', 'CAGCA', 'TGTAG', 'TGGAA', 'ACTAT', 'TACAG', 'TTAAG', 'TTGCG']
		#umi_list = ['AAAAA', 'TTTTT', 'CCCCC', 'GGGGG']
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)
		file_fq1 = gzip.open(self.fq1,'rb')
		file_fq2 = gzip.open(self.fq2,'rb')
		file_out_fq1 = open('%s/%s_1.splitUMI.fq'%(self.outdir,self.outprefix),'w')
		file_out_fq2 = open('%s/%s_2.splitUMI.fq'%(self.outdir,self.outprefix),'w')
		file_log = open('%s/%s.umi.log'%(self.outdir,self.outprefix),'w')
		umi_dict = {}
		for i in umi_list:
			umi_dict[i] = 0
		total_n = 0
		umi_n = 0
		NA_umi_n = 0
		umi_start = 0
		umi_end = 5
		for bline in file_fq1:
			total_n += 1
			line = bline.decode('ascii')
			seq = file_fq1.readline().decode('ascii')
			sig = file_fq1.readline().decode('ascii')
			qua = file_fq1.readline().decode('ascii')
			read2_id = file_fq2.readline().decode('ascii')
			read2_seq = file_fq2.readline().decode('ascii')
			read2_sig = file_fq2.readline().decode('ascii')
			read2_qua = file_fq2.readline().decode('ascii')
			umi1 = seq[umi_start:umi_end]
			umi2 = read2_seq[umi_start:umi_end]
			umi1_check = self.check_umi(umi_list,umi1,self.err)
			umi2_check = self.check_umi(umi_list,umi2,self.err)
			if line.split('/')[0] != read2_id.split('/')[0]:
				print('WARNING: The input fastq files are broken.Invalid fastq pair "%s | %s",start gzip ...'%(line.strip(),read2_id.strip()))
				break
			if umi1_check == 'NA' or umi2_check == 'NA':
				NA_umi_n += 1
				continue
			umi_n += 1
			umi_dict[umi1_check] += 1
			umi_dict[umi2_check] += 1
			umi_total = umi1_check + umi2_check
			read1_new_id = line.split('/')[0]+'_UMI_'+umi_total+'/1'
			read1_new_seq = seq[umi_end:]
			read1_new_qua = qua[umi_end:]
			read2_new_id = read2_id.split('/')[0]+'_UMI_'+umi_total+'/2'
			read2_new_seq = read2_seq[umi_end:]
			read2_new_qua = read2_qua[umi_end:]
			file_out_fq1.write(read1_new_id+'\n'+read1_new_seq+sig+read1_new_qua)
			file_out_fq2.write(read2_new_id+'\n'+read2_new_seq+sig+read2_new_qua)
		umi_rate = str('%.2f'%(100*umi_n/total_n))+'%'
		file_log.write('UMI_split_rate:\t'+umi_rate+'\n\n')
		file_log.close()
		#for key, value in umi_dict.items():
		#	rate = str('%.2f'%(100*int(value)/(total_n*2)))+'%'
		#	file_log.write(key+'\t'+str(value)+'\n')

		file_fq1.close()
		file_fq2.close()
		file_out_fq1.close()
		file_out_fq2.close()
		if os.path.exists('%s/%s_1.splitUMI.fq.gz'%(self.outdir,self.outprefix)):
			os.system('rm %s/%s_1.splitUMI.fq.gz'%(self.outdir,self.outprefix))
		if os.path.exists('%s/%s_2.splitUMI.fq.gz'%(self.outdir,self.outprefix)):
			os.system('rm %s/%s_2.splitUMI.fq.gz'%(self.outdir,self.outprefix))

		p = Pool(2)
		outfq1 = '%s/%s_1.splitUMI.fq'%(self.outdir,self.outprefix)
		outfq2 = '%s/%s_2.splitUMI.fq'%(self.outdir,self.outprefix)
		log = '%s/%s.umi.log'%(self.outdir,self.outprefix)
		p.starmap_async(self.gzip_fq,[(outfq1,log)])
		p.starmap_async(self.gzip_fq,[(outfq2,log)])
		p.close()
		p.join()

		with open(log,'a') as fl:
			fl.write('[ '+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+' ] '+'All gzip done.\n')
		return

class BamConsensusAnalysis(object):
	def __init__(self):
		self.inbam = optlist.opt_i
		self.outprefix = optlist.opt_o
		self.outdir = optlist.opt_d
		self.bed = optlist.opt_b
		self.familysize = optlist.opt_s
		self.mapq = optlist.opt_m
		self.ref = optlist.opt_F
		self.umi_type = optlist.opt_u
		self.base_err = optlist.opt_B
		if optlist.opt_r:
			self.rm = optlist.opt_r
		else:
			self.rm = 1
		if optlist.opt_X:
			self.java_xmx = optlist.opt_X
		else:
			self.java_xmx = '10G'
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)

	def GenerateTargetBam(self):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Generate target bam ')
		cmd = '%(samtools)s view -L %(bed)s -b %(inbam)s > %(outdir)s/%(outprefix)s.target.bam'\
			%{'samtools':samtools,'bed':self.bed,'inbam':self.inbam,'outdir':self.outdir,'outprefix':self.outprefix}
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Generate target bam finished ')
		return

	def RevertSam(self):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Revert bam ')
		cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(picard)s RevertSam SANITIZE=true I=%(outdir)s/%(outprefix)s.target.bam O=%(outdir)s/%(outprefix)s.target.revert.bam RESTORE_ORIGINAL_QUALITIES=false REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false MAX_DISCARD_FRACTION=1'\
			%{'java':java,'outdir':self.outdir,'outprefix':self.outprefix,'picard':picard,'xmx':self.java_xmx}
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] revert bam finished ')
		return

	def SetRX(self, inbam, outbam):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Set bam RX tags ')
		file_inbam = pysam.AlignmentFile(inbam,'rb',check_sq=False)
		header = file_inbam.header
		file_outbam = pysam.AlignmentFile(outbam,'wb',header=header)
		for r in file_inbam:
			readid = r.query_name
			UMI = readid.split('_')[-1]
			UMI_tag = UMI[0:int(len(UMI)/2)]+'-'+UMI[int(len(UMI)/2):]
			r.set_tag('RX',UMI_tag,value_type='Z')
			file_outbam.write(r)
		file_inbam.close()
		file_outbam.close()
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Set bam RX tags finished ')
		return

	def SortBam(self, inbam, outbam):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Sort bam ')
		cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s SortBam --input=%(inbam)s --output=%(outbam)s --sort-order=Queryname'\
			%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'xmx':self.java_xmx}
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Sort bam finished ')
		return

	def SetMate(self, inbam, outbam):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Set bam mate tags ')
		cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s SetMateInformation --input=%(inbam)s --output=%(outbam)s --ref=%(reference)s'\
			%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'reference':self.ref,'xmx':self.java_xmx}
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Set bam mate tags finished ')
		return

	def GroupUMI(self, inbam, outbam, mapq, umi_type):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Group UMI ')
		if umi_type == 'Duplex':
			cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s GroupReadsByUmi --input=%(inbam)s --output=%(outbam)s --min-map-q=%(mapq)s --edits=0 --raw-tag=RX --family-size-histogram=%(outdir)s/%(outprefix)s.family_size.xls --strategy=paired'\
				%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'mapq':self.mapq,'outdir':self.outdir,'outprefix':self.outprefix,'xmx':self.java_xmx}
		elif umi_type == 'Single':
			cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s GroupReadsByUmi --input=%(inbam)s --output=%(outbam)s --min-map-q=%(mapq)s --edits=0 --raw-tag=RX --family-size-histogram=%(outdir)s/%(outprefix)s.family_size.xls --strategy=adjacency'\
				%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'mapq':self.mapq,'outdir':self.outdir,'outprefix':self.outprefix,'xmx':self.java_xmx}
		else:
			print('ERROR: Invalid argument of umitype.Input Duplex/Single.')
			usage()
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Group UMI finished ')
		return

	def CallConsensus(self, inbam, outbam, familysize, umi_type):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Call consensus reads ')
		if umi_type == 'Duplex':
			cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s CallDuplexConsensusReads --min-reads=%(familysize)s --min-input-base-quality=0 --input=%(inbam)s --output=%(outbam)s'\
				%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'familysize':self.familysize,'xmx':self.java_xmx}
		elif umi_type == 'Single':
			cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s CallMolecularConsensusReads --min-reads=%(familysize)s --min-input-base-quality=0 --input=%(inbam)s --output=%(outbam)s'\
				%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'familysize':self.familysize,'xmx':self.java_xmx}
		else:
			print('ERROR: Invalid argument of umitype.Input Duplex/Single.')
			usage()
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Call consensus reads finished ')
		return

	def ReMap(self, inbam, outbam):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Remap bam ')
		cmd = '%(samtools)s fastq %(inbam)s | %(bwa)s mem -M -t 5 -p %(reference)s /dev/stdin | %(samtools)s view -F 256 -b - > %(outbam)s'\
			%{'samtools':samtools,'inbam':inbam,'outbam':outbam,'bwa':bwa,'reference':self.ref}
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Remap bam finished ')
		return

	def MergeBam(self, unmapbam, alignbam, outbam):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Merge bam files ')
		cmd = '%(java)s -Xmx%(xmx)s -jar -Djava.io.tmpdir=`pwd`/tmp %(picard)s MergeBamAlignment R=%(reference)s UNMAPPED_BAM=%(unmapbam)s ALIGNED_BAM=%(alignbam)s O=%(outbam)s CREATE_INDEX=true MAX_GAPS=-1 ALIGNER_PROPER_PAIR_FLAGS=true VALIDATION_STRINGENCY=SILENT SO=coordinate ATTRIBUTES_TO_RETAIN=XS'\
			%{'java':java,'picard':picard,'reference':self.ref,'unmapbam':unmapbam,'alignbam':alignbam,'outbam':outbam,'xmx':self.java_xmx}
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Merge bam files finished ')
		return

	def FilterConsensus(self, inbam, outbam, familysize, base_err, umi_type):
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] UMI correction ')
		if umi_type == 'Duplex':
			cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s FilterConsensusReads --input=%(inbam)s --output=%(outbam)s --ref=%(reference)s --min-reads=%(familysize)s --max-read-error-rate=0.025 --max-base-error-rate=%(base_err)s --min-base-quality=0 --max-no-call-fraction=0.20 --reverse-per-base-tags=True --require-single-strand-agreement=True'\
				%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'familysize':self.familysize,'reference':self.ref,'base_err':base_err,'xmx':self.java_xmx}
		elif umi_type == 'Single':
			cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s FilterConsensusReads --input=%(inbam)s --output=%(outbam)s --ref=%(reference)s --min-reads=%(familysize)s --max-read-error-rate=0.025 --max-base-error-rate=%(base_err)s --min-base-quality=0 --max-no-call-fraction=0.20 --reverse-per-base-tags=True --require-single-strand-agreement=false'\
				%{'java':java,'fgbio':fgbio,'inbam':inbam,'outbam':outbam,'familysize':self.familysize,'reference':self.ref,'base_err':base_err,'xmx':self.java_xmx}
		else:
			print('ERROR: Invalid argument of umitype.Input Duplex/Single.')
			usage()
		check_call(cmd,shell=True)
		print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] UMI correction finished ')
		return

	def main(self):
		target_bam = '%s/%s.target.bam'%(self.outdir,self.outprefix)
		revert_bam = '%s/%s.target.revert.bam'%(self.outdir,self.outprefix)
		RX_bam = '%s/%s.target.revert.RX.bam'%(self.outdir,self.outprefix)
		sort_bam = '%s/%s.target.revert.RX.sort.bam'%(self.outdir,self.outprefix)
		mate_bam = '%s/%s.target.revert.RX.sort.mate.bam'%(self.outdir,self.outprefix)
		group_bam = '%s/%s.target.revert.RX.sort.mate.group.bam'%(self.outdir,self.outprefix)
		consensus_ubam = '%s/%s.target.revert.RX.sort.mate.group.consensus.ubam'%(self.outdir,self.outprefix)
		consensus_bam = '%s/%s.target.revert.RX.sort.mate.group.consensus.bam'%(self.outdir,self.outprefix)
		merge_bam = '%s/%s.target.revert.RX.sort.mate.group.consensus.merge.bam'%(self.outdir,self.outprefix)
		merge_bam_bai = '%s/%s.target.revert.RX.sort.mate.group.consensus.merge.bai'%(self.outdir,self.outprefix)
		#filter_bam = '%s/%s.target.revert.RX.sort.mate.group.consensus.merge.filter.bam'%(self.outdir,self.outprefix)
		filter_bam = '%s/%s.target.consensus.bam'%(self.outdir,self.outprefix)
		self.GenerateTargetBam()
		self.RevertSam()
		self.SetRX(revert_bam, RX_bam)
		self.SortBam(RX_bam, sort_bam)
		self.SetMate(sort_bam,mate_bam)
		self.GroupUMI(mate_bam, group_bam, self.mapq, self.umi_type)
		self.CallConsensus(group_bam, consensus_ubam, self.familysize, self.umi_type)
		self.ReMap(consensus_ubam, consensus_bam)
		self.MergeBam(consensus_ubam, consensus_bam, merge_bam)
		self.FilterConsensus(merge_bam, filter_bam, self.familysize, self.base_err, self.umi_type)
		if int(self.rm) == 1:
			cmd = 'rm %s %s %s %s %s %s %s %s %s %s'%(target_bam,revert_bam,RX_bam,sort_bam,mate_bam,group_bam,consensus_ubam,consensus_bam,merge_bam,merge_bam_bai)
			check_call(cmd,shell=True)
		return

def usage():
	"""
UMI analysis package.
---------------------

Usage:
	UAP <command> [arguments]

Command list:
	AnnoFastqWithUMI <Annotates fastq file with UMI informations>
		Arguments:
		-f <FILE>   Input fastq1 file
		-r <FILE>   Input fastq2 file
		-e <INT>	Maximum mismatch bases in UMIs[Default:0]
		-o <STR>	Output prefix
		-d <STR>	Output directory

	BamConsensusAnalysis <Annotates existing BAM files and filter consensus reads>
		Arguments:
		-i <FILE>   Input bam file
		-o <STR>	Output prefix
		-d <STR>	Output directory
		-b <FILE>   Target bed file
		-m <INT>	Minimum mapping quality
		-s <INT>	Family size threshold
		-u <STR>	UMI type [Duplex/Single]
		-B <FLOAT>  Maximum base error rate in UMI family
		-F <STR>	Reference path
		-r <INT>	[0/1] 1 means remove all Temp files. Default:1

	ErrorRateStats <Calculates the error rate by read position on coordinate sorted mapped BAMs>
		Arguments:
		-i <FILE>   Input bam file
		-v <VCF>	Optional VCF file of variant sites to ignore
		-F <STR>	Reference path
		-m <INT>	Minimum mapping quality
		-o <STR>	Output prefix
	"""
	print(usage.__doc__)
	sys.exit(1)
	return

def ErrorRateStats(inbam, outprefix, ref, vcf, mapq, xmx):
	if vcf:
		cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s ErrorRateByReadPosition --input=%(inbam)s --output=%(outprefix)s --ref=%(ref)s --variants=%(vcf)s --min-mapping-quality=%(mapq)s'\
			%{'java':java,'fgbio':fgbio,'inbam':inbam,'outprefix':outprefix,'ref':ref,'vcf':vcf,'mapq':mapq,'xmx':xmx}
	else:
		cmd = '%(java)s -Xmx%(xmx)s -Djava.io.tmpdir=`pwd`/tmp -jar %(fgbio)s ErrorRateByReadPosition --input=%(inbam)s --output=%(outprefix)s --ref=%(ref)s --min-mapping-quality=%(mapq)s'\
			%{'java':java,'fgbio':fgbio,'inbam':inbam,'outprefix':outprefix,'ref':ref,'mapq':mapq,'xmx':xmx}
	check_call(cmd,shell=True)
	return

if __name__ == '__main__':
	if len(sys.argv) < 2:
		usage()
	parser = OptionParser()
	parser.add_option('-f', dest = 'opt_f', type = 'string')
	parser.add_option('-r', dest = 'opt_r', type = 'string')
	parser.add_option('-e', dest = 'opt_e', type = 'int')
	parser.add_option('-o', dest = 'opt_o', type = 'string')
	parser.add_option('-d', dest = 'opt_d', type = 'string')
	parser.add_option('-i', dest = 'opt_i', type = 'string')
	parser.add_option('-b', dest = 'opt_b', type = 'string')
	parser.add_option('-m', dest = 'opt_m', type = 'string')
	parser.add_option('-s', dest = 'opt_s', type = 'string')
	parser.add_option('-u', dest = 'opt_u', type = 'string')
	parser.add_option('-B', dest = 'opt_B', type = 'string')
	parser.add_option('-F', dest = 'opt_F', type = 'string')
	parser.add_option('-v', dest = 'opt_v', type = 'string')
	parser.add_option('-X', dest = 'opt_X', type = 'string')
	optlist, args = parser.parse_args()

	#binpath = os.path.dirname(sys.argv[0])
	#rootpath = os.path.dirname(binpath)
	UAP_HOME = os.environ.get('UAP_HOME')
	samtools = '%s/tools/samtools'%(UAP_HOME)
	bwa = '%s/tools/bwa'%(UAP_HOME)
	#java = '%s/tools/jre1.8.0_101/bin/java'%(rootpath)
	java = os.environ.get('JAVA_HOME')+'/bin/java'
	picard = '%s/tools/picard.jar'%(UAP_HOME)
	fgbio = '%s/tools/fgbio-1.0.0.jar'%(UAP_HOME)
	if 'help' in args:
		usage()
	elif 'AnnoFastqWithUMI' in args:
		annofastq = AnnoFastqWithUMI()
		annofastq.annofq()
	elif 'BamConsensusAnalysis' in args:
		bamconsensus = BamConsensusAnalysis()
		bamconsensus.main()
	elif 'ErrorRateStats' in args:
		if not optlist.opt_m:
			mapq = 20
		else:
			mapq = optlist.opt_m
		if optlist.opt_X:
			ErrorRateStats(optlist.opt_i,optlist.opt_o,optlist.opt_F,optlist.opt_v,mapq,optlist.opt_X)
		else:
			ErrorRateStats(optlist.opt_i,optlist.opt_o,optlist.opt_F,optlist.opt_v,mapq,'10G')
	else:
		usage()

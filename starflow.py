#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
import os, sys
import argparse
import hashlib
import glob
import uuid 
import re
import datetime


def get_parser():
	"""
	Parse arguments
	@return: arguments list
	@rtype: parser object
	"""

	parser = argparse.ArgumentParser(description='LRUE staph prestation')

	parser.add_argument('-i', action="store", dest='archive_path',
						type=str, required=True, help='archive path (REQUIRED)') 

	parser.add_argument('-p', action="store", dest='archive_password',
						type=str, required=True, help='archive password(REQUIRED)')

	parser.add_argument('-t', action="store", dest='resulTab',
						type=str, required=True, help='result tsv file (REQUIRED)')
						
	parser.add_argument('-RefCall', action="store", dest='RefCall',
						type=str, required=True, help='reference used for variant calling (REQUIRED)')

	parser.add_argument('-RefDir', action="store", dest='RefDir',
						type=str, required=True, help='references fasta files directory (REQUIRED)')
						
	parser.add_argument('-q', action="store", dest='NAuRA_queries',
						type=str, required=True, help='NAuRA query text file (REQUIRED)')					
						
	parser.add_argument('--ad', action="store", dest='ADAPTATERS',
						type=str, default="/global/bio/data/all_known_iontorrent-Illumina.fa", help='Sequencing adaptaters file (default:/global/bio/data/all_known_iontorrent-Illumina.fa)')

	parser.add_argument('--s', action="store", dest='xlsx2csv',
						type=str, default="/global/bio/bin/xlsx2csv", help='xlsx2csv path (default:/global/bio/bin/xlsx2csv)')

	parser.add_argument('--T', action="store", dest='nbThreads',
						type=str, default='8', help='number of threads (default:8)')	
						
	parser.add_argument('--SNVthreshold', action="store", dest='NumContamSNVs_threshold',
						type=int, default=2, help='NumContamSNVs threshold for confindr (default:2)')
						
	parser.add_argument('--minCov', action="store", dest='minCov',
						type=int, default=30, help='minimum coverage required for workflow (default:30)')
						
	parser.add_argument('--targetCov', action="store", dest='targetCov',
						type=int, default=100, help='minimum coverage required for normalization, 0 for no normalization (default:100)')
						
	parser.add_argument('--l', action="store", dest='minContigLen',
						type=int, default=200, help='minimum contig length (default:200)')
						
	parser.add_argument('--p', action="store", dest='minPhredScore',
						type=int, default=30, help='minimum phred score (default:30)')
						
	parser.add_argument('--a', action="store", dest='min_mapped_scaffold_percent', 
						type=int, default=80, help='minimum percent of total scaffold length \
						mapped on reference, used to detect contaminants (default:80)')
						
	parser.add_argument('--d', action="store", dest='max_len_diff', 
						type=int, default=30, help='maximum percent of length difference \
						between assembly and reference (default:20)')
						
	parser.add_argument('--Ksize', action="store", dest='Ksize', 
						type=str, default='15', help='Kmer size for mash (default:15)')	

	parser.add_argument('--Ssize', action="store", dest='Ssize', 
						type=str, default='1000', help='Sketch size for mash (default:15)')								

	return parser
	

class sample(object) :

	def __init__(self, SampleID, fastq1, fastq2, md5, csvFile, context, countrie, year, sequencing_Techno):

		self.id = SampleID	
		self.fastq1 = fastq1
		self.fastq2 = fastq2
		self.md5 = md5
		self.metadata = csvFile
		self.context = context
		self.countrie = countrie
		self.year = year
		self.sequencing_Techno = sequencing_Techno
		self.NumContamSNVs = "0"
		self.ContamStatus = ""
		self.PercentContam = "0"
		self.scaffoldingRef = ""
		
	def printItself(self):
		print("###################################")
		print("SampleID : " + self.id)
		print("fastq1 : " + self.fastq1)
		print("fastq2 : " + self.fastq2)
		print("md5 : " + self.md5)
		print("metadata : " + self.metadata)
		print("context : " + self.context)
		print("countrie : " + self.countrie)
		print("year : " + self.year)
		print("sequencing_Techno : " + self.sequencing_Techno)
		print("###################################")
		
	def writeItself(self,file):
		f=open(file,'a')
		f.write("###################################\n")
		f.write("SampleID : " + self.id + "\n")
		f.write("fastq1 : " + self.fastq1 + "\n")
		f.write("fastq2 : " + self.fastq2 + "\n")
		f.write("md5 : " + self.md5 + "\n")
		f.write("metadata : " + self.metadata + "\n")
		f.write("context : " + self.context + "\n")
		f.write("countrie : " + self.countrie + "\n")
		f.write("year : " + self.year + "\n")
		f.write("sequencing_Techno : " + self.sequencing_Techno + "\n")
		f.write("###################################" + "\n\n")
		f.close()
		
	def setContam(self,confindr_report,NumContamSNVs_threshold):
		f = open(confindr_report,'r')
		line = f.readlines()[1]
		f.close()
		self.NumContamSNVs = line.split(',')[2]
		if int(self.NumContamSNVs) <= NumContamSNVs_threshold :
			self.ContamStatus = False
		else:
			self.ContamStatus = True
		self.PercentContam = line.split(',')[4]
		
	def setRefScaffolding(self, ref_scaffolding):
		self.scaffoldingRef = ref_scaffolding

	def setArtworkParam(self,varcallRef,SNVthreshold,minCov,targetCov,minContigLen,minPhredScore,max_len_diff,Ksize,Ssize,min_mapped_scaffold_percent):
	
		self.varcallRef = varcallRef
		self.SNVthreshold = SNVthreshold
		self.minCov = minCov
		self.targetCov = targetCov
		self.minContigLen = minContigLen
		self.minPhredScore = minPhredScore
		self.max_len_diff = max_len_diff
		self.Ksize = Ksize
		self.Ssize = Ssize
		self.min_mapped_scaffold_percent = min_mapped_scaffold_percent
		
	def writeArtworkParam(self,file):
		f=open(file,'a')
		f.write("###################################\n")
		f.write("Variant calling reference : " + self.varcallRef + "\n")
		f.write("scaffolding reference : " + self.scaffoldingRef + "\n")
		f.write("SNV threshold : " + str(self.SNVthreshold) + "\n")
		f.write("min coverage : " + str(self.minCov) + "\n")
		f.write("target coverage (normalization) : " + str(self.targetCov) + "\n")
		f.write("min contig length : " + str(self.minContigLen) + "\n")
		f.write("min PhredScore : " + str(self.minPhredScore) + "\n")
		f.write("max length difference between assembly and scaffolding reference : " + str(self.max_len_diff) + "\n")
		f.write("min % identity between assembly and scaffolding reference : " + str(self.min_mapped_scaffold_percent) + "\n")
		f.write("Kmer size for mash : " + str(self.Ksize) + "\n")
		f.write("Sketch size for mash : " + str(self.Ssize) + "\n")
		f.write("###################################" + "\n\n")
		f.close()

	def compressReads(self):
		os.system("gzip " + self.fastq1)
		os.system("gzip " + self.fastq2)
		self.fastq1 = self.fastq1 + ".gz"
		self.fastq2 = self.fastq2 + ".gz"

	def setArtworkResult(self, result_folder):
		
		self.fastq1 = result_folder + os.path.basename(self.fastq1)
		self.fastq2 = result_folder + os.path.basename(self.fastq2)
		self.contig = result_folder + self.id + "_contig.fasta"
		self.assembly =  result_folder + self.id + "_assembly.fasta"
		self.vcf =  result_folder + self.id + ".vcf"
		self.mash =  result_folder + self.id + ".msh"
		self.Artwork_report =  result_folder + self.id + "_report.txt"
		self.Artwork_log =  result_folder + self.id + ".log"
		self.gbk =  result_folder + self.id + ".gbk"
		self.gff =  result_folder + self.id + ".gff"
		self.R1fastqc_data = result_folder + self.id + "_R1_fastqc/fastqc_data.txt"
		self.R1fastqc_per_base_quality = result_folder + self.id + "_R1_fastqc/Images/per_base_quality.png"
		self.R1fastqc_per_sequence_gc_content = result_folder + self.id + "_R1_fastqc/Images/per_sequence_gc_content.png"
		self.R1fastqc_per_sequence_quality = result_folder + self.id + "_R1_fastqc/Images/per_sequence_quality.png"
		self.R2fastqc_data = result_folder + self.id + "_R2_fastqc/fastqc_data.txt"
		self.R2fastqc_per_base_quality = result_folder + self.id + "_R2_fastqc/Images/per_base_quality.png"
		self.R2fastqc_per_sequence_gc_content = result_folder + self.id + "_R2_fastqc/Images/per_sequence_gc_content.png"
		self.R2fastqc_per_sequence_quality = result_folder + self.id + "_R2_fastqc/Images/per_sequence_quality.png"
		self.quast_report = result_folder + "QUAST/report.tsv"
		self.setSequenceType()
		self.setQuastStats()
		
	def setNauraResult(self,dico, query_list):
	
		self.naura_matrix = os.path.dirname(os.path.realpath(self.gbk)) + '/' + self.id + "_NAuRA.tsv"
		self.naura_newAlleleDico = dico
		self.naura_paramFile = query_list

	def setSequenceType(self):
	
		f = open(self.Artwork_report,'r')
		last_line = f.readlines()[-1]
		self.ST = last_line.split("'")[3]
	
	def setQuastStats(self):
		
		f = open(self.quast_report,'r')
		lines = f.readlines()
		f.close()
		
		nbContig = lines[1].rstrip().split('\t')[1]
		N50 = lines[19].rstrip().split('\t')[1]
		genomeSize = lines[15].rstrip().split('\t')[1]
		
		self.nbContig = nbContig
		self.N50 = N50
		self.genomeSize = genomeSize


def extract(archive_path, archive_password):

	folder_path = os.path.basename(archive_path).split('.')[0]
	commande = "unzip -P " + archive_password + " " + archive_path + " -d " + folder_path
	os.system(commande)

	return folder_path
	
	
def xlsx_to_csv(xlsx_file, log, xlsx2csv):	

	try :
		open(xlsx_file,'r')
	except :
		logFile = open(log,'a')
		logFile.write("ERROR :: " + xlsx_file + " not found !\n\n")
		logFile.close()
		sys.exit(1)

	xlsx_filename = os.path.basename(xlsx_file)
	xlsx_dir = os.path.dirname(os.path.realpath(xlsx_file))
	xlsx_filename_without_suffix = xlsx_filename.split('.')[0]

	csv_path = xlsx_dir + '/' + xlsx_filename_without_suffix + ".csv"
	commande = xlsx2csv + " " + xlsx_file + " " + csv_path
	os.system(commande)

	return csv_path
	
	
def make_md5(archive_path):

	BLOCKSIZE = 65536
	hasher = hashlib.md5()

	with open(archive_path, 'rb') as afile:
		buf = afile.read(BLOCKSIZE)
		while len(buf) > 0:
			hasher.update(buf)
			buf = afile.read(BLOCKSIZE)

	md5 = hasher.hexdigest()
	
	return md5
	
	
def grep_md5_in_result_file(md5, resulTab):

	resulTab_file = open(resulTab, 'r')
	lines = resulTab_file.readlines()
	resulTab_file.close()
	
	for line in lines :
		if md5 in line :
			return True
			
	return False


def read_csv_metadata(csv_file, md5, folder, logFile):

	csv = open(csv_file,'r')
	lines = csv.readlines()
	csv.close()
	
	list_result = []
	header = True
	
	for line in lines :
		if header :
			header = False
			continue
		line = line.rstrip().split(',')
		if line[0] == '':
			continue
			
		context = line[0]
		sampleID = line[1]
		countrie = line[2]
		year = line[3]
		sequencing_Techno = line[4]
		fastq1 = line[5]
		fastq2 = line[6]
		
		check_files_in_folder(folder, [fastq1,fastq2], logFile)
		fastq1 = os.path.dirname(os.path.realpath(csv_file)) + '/' + fastq1
		fastq2 = os.path.dirname(os.path.realpath(csv_file)) + '/' + fastq2
		new_sample = sample(sampleID, fastq1, fastq2, md5, csv_file, context, countrie, year, sequencing_Techno)
		list_result.append(new_sample)

	return list_result


def check_files_in_folder(folder, files_list, logFile):

	list_file_in_folder = glob.glob(folder + "/*")
	list_file_in_folder_without_path = []
	
	for element in list_file_in_folder :
		list_file_in_folder_without_path.append(os.path.basename(element))
	
	flag = True
	log = open(logFile,'a')
	for element in files_list :
		if element not in list_file_in_folder_without_path :
			flag = False
			log.write("ERROR :: file " + element + " not found in " + folder + "\n")
	
	log.close()
	if flag == False :
		sys.exit(1)
	
	return True


def contamfinder(sampleObj, nbThreads):

	reads_folder = os.path.dirname(os.path.realpath(sampleObj.fastq1))
	confindr_output_folder = os.path.dirname(os.path.realpath(sampleObj.fastq1)) + "/confindr"
	
	confindr_command = "confindr.py -i " + reads_folder + ' -t ' + nbThreads + ' -o ' + confindr_output_folder
	command = '/bin/bash -c "source /global/conda/bin/activate confindr; ' + confindr_command + '"'
	os.system(command)
	
	report_file = confindr_output_folder + "/confindr_report.csv"
	
	return report_file
	
	
def scafolding_ref_finder(sampleObj, RefDir, nbThreads):	

	list_file_in_folder = glob.glob(RefDir + "/*.fa*")
	tmpFile = str(uuid.uuid4())
	fastqTmp = str(uuid.uuid4())
	fastqMsh = str(uuid.uuid4())
	
	cat_command = 'cat ' + sampleObj.fastq1 + ' ' + sampleObj.fastq2 + ' > ' + fastqTmp
	os.system(cat_command)
	
	mash_command = "mash sketch -o " + tmpFile + ' ' + " ".join(list_file_in_folder)
	command = '/bin/bash -c "source /global/conda/bin/activate artwork; ' + mash_command + '"'
	os.system(command)
	
	mash_command = "mash sketch -o " + fastqMsh + ' ' + fastqTmp
	command = '/bin/bash -c "source /global/conda/bin/activate artwork; ' + mash_command + '"'
	os.system(command)
	os.system("rm " + fastqTmp)
	
	mash_dist_output = os.path.dirname(os.path.realpath(sampleObj.fastq1)) + "/mash_dist.tsv"
	mash_command = "mash dist -p " + str(nbThreads) + ' ' + fastqMsh + ".msh " + tmpFile + ".msh > " + mash_dist_output
	command = '/bin/bash -c "source /global/conda/bin/activate artwork; ' + mash_command + '"'
	os.system(command)
	os.system("rm " + tmpFile + ".msh " + fastqMsh + ".msh")

	mash_resultFile = open(mash_dist_output,'r')
	lines = mash_resultFile.readlines()
	mash_resultFile.close()

	genome_ref = ""
	minDist = 1.0

	for line in lines :
		
		line = line.rstrip().split('\t')
		genome = line[1]
		dist = float(line[2])

		if dist < minDist :
			minDist = dist
			genome_ref = genome
			
	cp_command = "cp " + genome_ref + " " + os.path.dirname(os.path.realpath(sampleObj.fastq1))
	os.system(cp_command)
	ref_scaffolding = os.path.dirname(os.path.realpath(sampleObj.fastq1)) + '/' + os.path.basename(genome_ref)
	
	os.system("rm " + mash_dist_output)
	sampleObj.setRefScaffolding(ref_scaffolding)
			

def ARTwork_run(sampleObj, nbThreads, adaptaters):

	workdir = os.path.dirname(os.path.realpath(sampleObj.fastq1))
	os.mkdir(workdir + "/artwork")
	workdir = workdir + "/artwork"
	
	sampleObj.compressReads()
	
	### c'est moche
	cp_command = "mv " + sampleObj.fastq1 + " " + os.path.dirname(os.path.realpath(sampleObj.fastq1)) + "/artwork/."
	os.system(cp_command)
	cp_command = "mv " + sampleObj.fastq2 + " " + os.path.dirname(os.path.realpath(sampleObj.fastq1)) + "/artwork/."
	os.system(cp_command)
	
	fastq1 = os.path.basename(sampleObj.fastq1)
	fastq2 = os.path.basename(sampleObj.fastq2)
	###
	
	artwork_command = "ARtWORK_LRUE -1 " + fastq1 + ' -2 ' + fastq2 + ' -r1 ' + sampleObj.varcallRef + ' -r2 ' + sampleObj.scaffoldingRef + ' -ad ' + adaptaters + ' -minCov ' + str(sampleObj.minCov) + ' -targetCov ' + str(sampleObj.targetCov) + ' -l ' + str(sampleObj.minContigLen) + ' -p ' + str(sampleObj.minPhredScore) + ' -a ' + str(sampleObj.min_mapped_scaffold_percent) + ' -Ksize ' + str(sampleObj.Ksize) + ' -Ssize ' + str(sampleObj.Ssize) + ' -d ' + str(sampleObj.max_len_diff) + ' -T '  +  str(nbThreads) + ' -w ' + workdir
	
	command = '/bin/bash -c "source /global/conda/bin/activate artwork; ' + artwork_command + '"'
	os.system(command)
	
	result_folder = os.path.dirname(os.path.realpath(sampleObj.fastq1)) + "/artwork/artwork_results_for_" + sampleObj.id + "/"
	sampleObj.setArtworkResult(result_folder)


def ARTwork_grep_errors(sampleObj,logFile):

	report = open(sampleObj.Artwork_report, 'r')
	lines = report.readlines()
	report.close()
	
	warning1 = "ERROR : assembly length too different from"
	warning2 = "ERROR : only "
	warning3 = "low coverage : "
	
	for line in lines :
		if (warning1 in line) or (warning2 in line) or (warning3 in line):
			log = open(logFile,'a')
			log.write(line)
			log.close()
			return False

	return True


def NAuRA_run(sampleObj, query_list, nbThreads):

	queries_before_NAuRA = NAuRA_reads_queries_fasta(query_list)

	gbk_folder = os.path.dirname(os.path.realpath(sampleObj.gbk))
	naura_command = "NAuRA -i " + gbk_folder + " -q " + query_list + " -T " + nbThreads + " --noDrift"
	command = '/bin/bash -c "source /global/conda/bin/activate naura; ' + naura_command + '"'
	os.system(command)
	
	queries_after_NAuRA = NAuRA_reads_queries_fasta(query_list)
	dico_diff = NAuRA_compare_dico(queries_before_NAuRA, queries_after_NAuRA)
	sampleObj.setNauraResult(dico_diff,query_list)
	
	mv_command = "mv matrix.tsv " + gbk_folder + '/' + sampleObj.id + "_NAuRA.tsv"
	os.system(mv_command)
	os.system("rm list.txt")

	return
	
	
def fasta_to_dico(fasta):

	dico_fasta = {}
	fasta_file = open(fasta,'r')
	lines = fasta_file.readlines()
	seqId = ""
	sequence = ""
	for line in lines :
		if len(line) == 0 :
			continue
		if line[0] == '>' :
			if seqId != "" :
				dico_fasta[seqId] = sequence
			seqId = line.rstrip()[1:]
			sequence = ""
		else :
			sequence = sequence + line.rstrip()
	dico_fasta[seqId] = sequence
			
	return dico_fasta

	
def NAuRA_reads_queries_fasta(query_list):

	dico_query = {}
	qFile = open(query_list,'r')
	lines = qFile.readlines()
	qFile.close()
	
	for line in lines :
		fasta_file = line.split('\t')[0]
		dico_query[NAuRA_query_name(fasta_file)] = fasta_to_dico(fasta_file)
		
	return dico_query


def NAuRA_query_name(fasta):

	f=open(fasta,'r')
	lines = f.readlines()
	f.close()
	for line in lines :
		if line[0] == '>' :
			return '_'.join(line.rstrip()[1:].split('_')[0:-1])


def NAuRA_combin_result(sampleObj_list, folder) :

	if len(sampleObj_list) == 1 :
		return sampleObj_list[0].naura_matrix
		
	matrix_combin = folder + "/NAuRA_result_combin.tsv"
	matrix_combin_file = open(matrix_combin,'w')
	firstObj = True
	
	for sample in sampleObj_list :
		f = open(sample.naura_matrix,'r')
		lines = f.readlines()
		f.close()
		if firstObj :
			firstObj = False
			matrix_combin_file.write(lines[0])
		if lines[1][-1]!='\n':
			lines[1] = lines[1] + '\n'
		matrix_combin_file.write(lines[1])
	
	matrix_combin_file.close()
	return matrix_combin


def NAuRA_compare_dico(dico_before,dico_after):

	dico_diff = {}

	for query in dico_after :
		for allele in dico_after[query]:
			if allele not in dico_before[query]:
				dico_diff[allele] = dico_after[query][allele]
				print(allele)
				break

	return dico_diff


def write_summary_results(resultFile, sampleObj):

	f = open(resultFile,'a')
	result_list = []
	result_list.append(datetime.datetime.today().strftime('%Y-%m-%d'))
	result_list.append(sampleObj.context)
	result_list.append(sampleObj.id)
	result_list.append(sampleObj.countrie)
	result_list.append(sampleObj.year)
	result_list.append(sampleObj.sequencing_Techno)
	result_list.append(sampleObj.genomeSize)
	result_list.append(sampleObj.N50)
	result_list.append(sampleObj.nbContig)
	result_list.append(sampleObj.ST)
	result_list.append(parse_NAuRA_matrix_for_resultFile(sampleObj.naura_matrix))
	new_allele_dico = sampleObj.naura_newAlleleDico
	if len(new_allele_dico.keys())==0 :
		result_list.append("NA")
	else :
		result_list.append('|'.join(list(new_allele_dico.keys())))
	result_list.append(sampleObj.md5)
	f.write('\t'.join(result_list)+'\n')
	f.close()
	
	return
	
	
def parse_NAuRA_matrix_for_resultFile(NAuRA_matrix):

	f = open(NAuRA_matrix,'r')
	lines = f.readlines()
	f.close()
	
	SEs_list = lines[0].rstrip()[1:].split('\t')
	NAuRA_string = ""
	
	line = lines[1].rstrip().split('\t')[1:]
	x = 0
	for element in line :
		element = element.replace(' ','')
		if element != '0':
			if NAuRA_string != "" :
				NAuRA_string = NAuRA_string + '|'
			NAuRA_string = NAuRA_string + SEs_list[x] + '[' + element + ']'
		x+=1 
		
	return NAuRA_string

	
def sampleObj_filter(sample_list,sample_to_reject):

	sampleID_to_reject = []
	for sample in sample_to_reject :
		sampleID_to_reject.append(sample.id)
		
	sampleList_toKeep = []
	for sample in sample_list :
		if sample.id not in sampleID_to_reject :
			sampleList_toKeep.append(sample)
	
	return sampleList_toKeep
		

def make_kmer_tree(sampleObj_list, RefDir, nbThreads, folder):

	tmpFile = str(uuid.uuid4())
	f = open(tmpFile,'w')
	list_file_in_folder = glob.glob(RefDir + "/*")
	for element in list_file_in_folder :
		f.write(element + '\n')
	for sample in sampleObj_list :
		f.write(sample.assembly + '\n')
	f.close()
	
	output_matrix = "kmer_matrix"
	output_tree = "kmer_tree"

	fastosh_command = "FasTosh -i " + tmpFile + " -o " + output_matrix + " --S -T " + nbThreads + " -e " + output_tree
	command = '/bin/bash -c "source /global/conda/bin/activate fastosh; ' + fastosh_command + '"'
	os.system(command)
	
	os.system("rm " + tmpFile)
	
	f = open(output_tree + ".nwk",'r')
	line = f.readlines()[0]
	f.close()
	line = line.replace("'",'').replace(".fasta",'').replace("Staphylococcus_aureus_",'')
	f = open(output_tree + ".nwk",'w')
	f.write(line)
	f.close()
	
	os.system("mv " + output_matrix + ".tsv " + output_tree + ".nwk " + folder)
	output_matrix = folder + '/' + output_matrix + ".tsv"
	output_tree = folder + '/' + output_tree + ".nwk"
	output_pdf = folder + "/kmer_tree.pdf"
	
	figtree_command = "figtree -graphic PDF " + output_tree + " " + output_pdf
	command = '/bin/bash -c "source /global/conda/bin/activate fastosh; ' + figtree_command + '"'
	os.system(command)

	return output_matrix,output_tree,output_pdf
	

#main function	
def main():	

	XLS_FILENAME = "NAuRA_online_request.xlsx"
	
	####### variable pour test ####### 
	folder = "Archive_test"
	md5 = "238f28bb88e1562a9208173cc674d20e"
	csv_input = "Archive_test/NAuRA_online_request.csv"
	confindr_report = "Archive_test/confindr/confindr_report.csv"
	####### 

	##########################################
	#			Initialisation				 #
	##########################################

	
	# Get arguments 
	parser=get_parser()
	
	# Print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
		
	Arguments=parser.parse_args()
	
	logFile = os.path.basename(Arguments.archive_path).split('.')[0] + ".log"
	
	### make md5 from archive	
	log = open(logFile, 'a')
	log.write("# Create md5 from archive\n")
	md5 = make_md5(Arguments.archive_path)
	log.write("--> " + md5 + "\n\n")
	log.close()
	
	### check if md5 already in result tab
	log = open(logFile, 'a')
	log.write("# Check if archive already been analyzed\n")
	checkpoint = grep_md5_in_result_file(md5, Arguments.resulTab)
	if checkpoint :
		log.write("ERROR :: " + md5 + " exist in file " + Arguments.resulTab + "\n\n")
		sys.exit(1)
	else :
		log.write("pass\n\n")
	log.close()
	
	### extract archive
	log = open(logFile, 'a')
	log.write("# Archive extraction\n")
	folder = extract(Arguments.archive_path, Arguments.archive_password)
	log.write("Ok\n\n")
	log.close()
	
	### convert xlsx file
	log = open(logFile, 'a')
	log.write("# Convert xlsx to csv\n")
	xlsx_file = folder + '/' + XLS_FILENAME
	csv_input = xlsx_to_csv(xlsx_file,logFile, Arguments.xlsx2csv)
	log.write("Ok\n\n")
	log.close()
		
	### read input tab file, check files and make sample object
	log = open(logFile, 'a')
	log.write("# Check files in folder\n")
	log.close()
	sample_list = read_csv_metadata(csv_input, md5, folder, logFile)
	log = open(logFile, 'a')
	log.write("Ok\n\n")
	log.close()
	
	sample_to_reject = []
	
	for sample in sample_list:
		sample.writeItself(logFile)
		
		### Detect contaminant
		log = open(logFile, 'a')
		log.write("# Check reads contam\n")
		confindr_report = contamfinder(sample, Arguments.nbThreads)
		sample.setContam(confindr_report,Arguments.NumContamSNVs_threshold)
		log.write("NumContamSNVs = " + sample.NumContamSNVs + "\n")
		log.write("PercentContam = " + sample.PercentContam + "\n")
		log.write("ContamStatus = " + str(sample.ContamStatus) + "\n")
		if sample.ContamStatus :
			sample_to_reject.append(sample)
			log.write("Error :: " + sample.id + " contaminated ! Analyses stopped for this sample.\n")
			continue
		log.close()
		
		### mash for closest ref
		log = open(logFile, 'a')
		log.write("\n# Search closest reference with mash\n")
		scafolding_ref_finder(sample, Arguments.RefDir, Arguments.nbThreads)
		log.write("--> " + os.path.basename(sample.scaffoldingRef) + "\n")
		log.close()
		#
		#sample.scaffoldingRef = "/global/scratch/a.felten/LRUE_Staph/prestation_scripts/Archive_test/Staphylococcus_aureus_FR_821779.fasta"
		
		### artwork
		log = open(logFile, 'a')
		log.write("\n# Run ARTwork\n")
		log.close()	
		sample.setArtworkParam(Arguments.RefCall,Arguments.NumContamSNVs_threshold,Arguments.minCov,Arguments.targetCov,Arguments.minContigLen,Arguments.minPhredScore,Arguments.max_len_diff,Arguments.Ksize,Arguments.Ssize, Arguments.min_mapped_scaffold_percent)
		sample.writeArtworkParam(logFile)
		ARTwork_run(sample, Arguments.nbThreads, Arguments.ADAPTATERS)
		pass_Artwork = ARTwork_grep_errors(sample,logFile)
		if not pass_Artwork :
			sample_to_reject.append(sample)
			continue
		log = open(logFile, 'a')
		log.write("Finish\n\n")
		log.close()
		
		#sample.gbk = "/global/scratch/a.felten/LRUE_Staph/prestation_scripts/Archive_test/artwork/artwork_results_for_05CEB51/05CEB51.gbk"
		
		### NAuRA
		log = open(logFile, 'a')
		log.write("\n# Run NAuRA\n")
		log.close()	
		NAuRA_run(sample, Arguments.NAuRA_queries, Arguments.nbThreads)
		#print(sample.naura_newAlleleDico)
		log = open(logFile, 'a')
		log.write("Finish\n\n")
		log.close()
		
		### write result in tab
		log = open(logFile, 'a')
		log.write("\n# write results in " + Arguments.resulTab + "\n\n")
		log.close()	
		write_summary_results(Arguments.resulTab, sample)
		
		sample.assembly = "Archive_test/artwork/artwork_results_for_05CEB51/05CEB51_assembly.fasta"
		
	### remove sample to reject
	if len(sample_to_reject)>0:
		sample_list = sampleObj_filter(sample_list,sample_to_reject)
	log = open(logFile, 'a')
	log.write("\n# Sample analyzed :\n")
	list_id = []
	for sample in sample_list :
		list_id.append(sample.id)
	log.write(" - ".join(list_id)+"\n\n")
	log.close()	
	
		
	### combin all NAuRA matrix
	log = open(logFile, 'a')
	log.write("\n# Combin all matrix from NAuRA output\n\n")
	log.close()	
	#NAuRA_combin_filePath = NAuRA_combin_result(sample_list, folder) 	
		
	### run fastosh
	log = open(logFile, 'a')
	log.write("\n# Compute clustering tree from all sample and reference genome\n")
	kmer_matrix, kmer_tree, kmer_pdf = make_kmer_tree(sample_list, Arguments.RefDir, Arguments.nbThreads, folder)
	log.write("Matrix file : " + kmer_matrix + "\n")
	log.write("Newick file : " + kmer_tree + "\n")
	log.write("PDF file : " + kmer_pdf + "\n\n")
	log.close()	
	
	# convert quast tsv to csv (only col 1+2)
	# tr '\t' ',' < report.tsv | sed 's/^# //g' | cut -d ',' -f 1,2 | sed 's/>=/$\\\geq$/g' | sed 's/%/\\\%/g' | sed 's/_/-/g' > report.csv
		
if __name__ == "__main__":
	main()
	
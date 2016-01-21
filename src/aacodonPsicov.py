#aacodonPsicov.py
# PSICOVcodon - Protein Sparse Inverse COVariance analysis program for nucleotide and amino acid sequences
# By Etai Jacob, September 2015 - Copyright (C) 2015 Etai Jacob, etai.jacob@gmail.com
# This code is licensed under the terms of GNU General Public License v2 or later

# This code was using in the following paper:

#Codon-level information improves predictions of inter-residue contacts in proteins by correlated mutation analysis

#Etai Jacob, Ron Unger, Amnon Horovitz
#Bar-Ilan University, Israel; Weizmann Institute of Science, Israel
#DOI: http://dx.doi.org/10.7554/eLife.08932
#Published September 15, 2015
#Cite as eLife 2015;4:e08932
#- See more at: http://elifesciences.org/content/4/e08932


#TODOS:
# 3. Implement psicov python warapper to aacodonPsicov together
# 4. Check if the parallel version can work with the compiler


#this script takes a fasta alignment file and an sth file and creates an aln file for psicov input based on indices exclusion of a reference seq in an sth file
#char excluded are '.' and lower case letters in the reference seq founded in the sth file
from Bio import AlignIO, SeqIO
import sys, re, numpy, getopt, os, sqlite3, csv

#Place the path to the PSICOV binary
codonPsicov_cmd = "/dir/leads_test/etai/SIMU/DB/PSICOV/PSICOV.V2.1b3/cmc -l -p -r 0.001" #-d 0.03" was -p -r 0.001, -l means no APC
aaPsicov_cmd = "/dir/leads_test/etai/SIMU/DB/PSICOV/PSICOV.V2.1b3/amc -l -p -r 0.001" # 0.03" (Recommended options for quicker results) was -p -r 0.001


def main(argv=None):
	if argv is None:
		argv = sys.argv

	AAfastaFile = None
	CfastaFile = None
	sthFile = None
	refName = None

	options, remainder = getopt.getopt(sys.argv[1:], 'a:c:s:r:', ['aa_msa=', 
                                                         	   	 'codon_msa=',
                                                         	   	 'sth_msa=',
                                                         	   	 'ref_name=', ])
	print 'OPTIONS   :', options
	print 'REMINDER  :', remainder
	for opt, arg in options:
		if opt in ('-a', '--aa_msa'):
			AAfastaFile = arg
		elif opt in ('-c', '--codon_msa'):
			CfastaFile = arg
		elif opt in ('-s', '--sth_msa'):
			sthFile = arg
		elif opt in ('-r', '--ref_name'):
			refName = arg

	if sthFile is None:
		print "Please enter an sthFile (-s or --sth_msa)"
		return 1
	# if refName is None:
	# 	print "Please enter a refName (-r or --ref_name)"
	# 	return 1
	if AAfastaFile is None and CfastaFile is None:
		print "Please enter a codonFile (-c or --codon_msa) or an AA file (-a or --aa_msa)"
		return 1
		

	refSeq = grepRefseq(sthFile, refName)
	idx2exclude = findSth2ExcludeOccurences(refSeq)
	print "Excluding %d positions from re seq name %s" % (len(idx2exclude), refName)
	if AAfastaFile is not None:
		#Write AA aln file:
		with open(AAfastaFile, "r") as f:
			with open(AAfastaFile + ".aln", "w") as fo:
				print "Writing AA fasta file."
				for record in SeqIO.parse(f, "fasta"):
					seqout = ''.join(numpy.delete(numpy.array(list(record.seq)), idx2exclude))
					fo.write(seqout + "\n")

		
	if CfastaFile is not None:
		#Write Codon aln file:
		with open(CfastaFile, "r") as f:
			with open(CfastaFile + ".aln", "w") as fo:
				print "Writing Codon fasta file."
				for record in SeqIO.parse(f, "fasta"):
					na = numpy.array([str(record.seq[(i*3):(i*3+3)]) for i in range(len(record.seq)/3)])
					seqout = ''.join(numpy.delete(na, idx2exclude))
					fo.write(seqout + "\n")
				
	#Run psicov for codon and aa:
	if AAfastaFile is not None:
		cmd = aaPsicov_cmd + " " + AAfastaFile + ".aln > " + AAfastaFile + ".aln.aapsicov.r.noAPC.out"
		print cmd
		os.system(cmd)

	if CfastaFile is not None:
		cmd = codonPsicov_cmd + " " + CfastaFile + ".aln > " + CfastaFile + ".aln.codonpsicov.r.noAPC.out"
		print cmd
		os.system(cmd)

	# dbconn = sqlite3.connect(":memory:")
	# dbc = dbconn.cursor()
	# acm = csv.reader(open(AAfastaFile + ".aln.aapsicov.out"), delimiter = " ")
	# ccm = csv.reader(open(CfastaFile + ".aln.codonpsicov.out"), delimiter = " ")
	# dbc.execute("create table acm (i int, j int, aascore float)")
	# dbc.execute("create table ccm (i int, j int, codonscore float)")
	# dbc.executemany("insert into acm(i, j, aascore) values (?, ?, ?)", acm)
	# dbc.executemany("insert into ccm(i, j, codonscore) values (?, ?, ?)", ccm)

	# dbc.execute("select * from acm limit 5")
	# print dbc.fetchone()
	# dbconn.close()


def findSth2ExcludeOccurences(s, ch = '.'):
    return [i for i, letter in enumerate(s) if letter == ch or letter.islower()]

def grepRefseq(sthFile, refName = None):
	refSeq = None

	with open(sthFile, "r") as f:
		for line in f: #SeqIO.parse(f, "stockholm"):
			line = line.rstrip()
			if len(line) < 3:
				continue
			if re.search("^#", line) is None:
				pl = line.split()
				if len(pl) < 2:
					#print "line is too short!!"
					continue
				if refName is None:
					refSeq = pl[1]
					print "Taking the first seq as a ref seq since no refName was given. \nseq: %s\nname: %s" % (refSeq, pl[0])
					break
				else:
					if pl[0] == refName:
						refSeq = pl[1]
						#print "Found! %s\n%s" % (refName, refSeq)
						break

	return refSeq


if __name__ == "__main__":
	sys.exit(main())

#/dir/leads_test/etai/SIMU/DB/DCA/cmds >run_commands.pl -port 19970 -cmds aacodonpsicov3.cmds -hosts 6@duper:1-4 -mail_to etai@cgen.com -cookie /dir/leads_test/etai/SIMU/DB/DCA/cmds/aacodonpsicov3.cmds.cookie -bad_file /dir/leads_test/etai/SIMU/DB/DCA/cmds/aacodonpsicov3.cmds.bad_file -status_file /dir/leads_test/etai/SIMU/DB/DCA/cmds/aacodonpsicov3.cmds.status_file > & ! aacodonpsicov3.cmds.runner &
#UnBatch_control -runner duper1 -rport 19971 -cookie aacodonpsicov3.cmds.cookie -hosts 8@duper:1-4

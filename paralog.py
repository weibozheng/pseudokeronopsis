import os
import sys,getopt

help_str='paralog.py -i <fasta file prot or nucl> -o <output folder> -d <prot or nucl>'
outfolder=''
fasta_file=''
data_type=''
try:
	opts,args=getopt.getopt(sys.argv[1:],"hi:d:o:",["help"])
except getopt.GetoptError:
	print (help_str)
	sys.exit(2)
for opt,value in opts:
	if opt in ("-h","--help"):
		print (help_str)
		sys.exit()
	if opt in ("-o"):
		outfolder=value
	if opt in ("-i"):
		fasta_file=value
	if opt in ("-d"):
		data_type=value
		if data_type !='prot' and data_type !='nucl':
			print("wrong data type, exit")
			sys.exit()
if not(outfolder and fasta_file and data_type):
	print (help_str)
	sys.exit()
else:
	os.system('mkdir '+ outfolder)

ol=os.listdir(outfolder)
if 'target.pdb' not in ol:
	os.system('makeblastdb -in '+fasta_file+' -dbtype '+data_type+' -out '+outfolder+'/target')
else:
	print('database already built, skip')
if 'target_self.tab' not in ol:
	if data_type=='prot':
		os.system("blastp -query "+fasta_file+' -db '+outfolder+'/target -outfmt 6 -num_threads 16 -max_target_seqs 10 -out '+outfolder+'/target_self.tab')
	else:
		os.system("blastn -query "+fasta_file+' -db '+outfolder+'/target -outfmt 6 -num_threads 16 -max_target_seqs 10 -out '+outfolder+'/target_self.tab')
else:
	print('self blast tabular file already built, skip')
inf=open(outfolder+'/target_self.tab')
dict_len=dict()
para_num=dict()
for line in inf:
	ll=line.split("\t")
	this_tag=ll[0]
	target_tag=ll[1]
	align_length=int(ll[3])
	if this_tag==target_tag:
		dict_len.update({this_tag:align_length})
inf=open(outfolder+'/target_self.tab')
for line in inf:
	ll=line.split("\t")
	this_tag=ll[0]
	target_tag=ll[1]
	iden=float(ll[2])
	align_length=int(ll[3])
	mismatch=int(ll[4])
	if (this_tag!=target_tag) and (iden>60 and align_length/dict_len[this_tag]>0.6):
		if this_tag not in para_num:
			para_num.update({this_tag:1})
		else:
			para_num[this_tag]+=1
para_count_file=outfolder+'/para_count.tab'
outf=open(para_count_file,"w")
for key in para_num:
	print(key+"\t"+str(para_num[key]),file=outf)
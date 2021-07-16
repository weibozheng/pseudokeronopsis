import MySQLdb
import re
import pandas
import rpy2
import math
import numpy
from Bio import Seq

#############parameter region##################
database_name='pse'
table_name1='gff_flava'
value_tag='sequence'

limit_rule1="where type='intron'"
coord_limit=""
out1=open(r"C:\Users\wbzhe\Pse\pflava_introns_weblog.txt","w")
###############################################
db=MySQLdb.connect("localhost", "root", "", database_name, charset='utf8' )
command_txt1="select "+value_tag+" from "+table_name1+" "+limit_rule1+";"
results1=pandas.read_sql(command_txt1,db)
for i in results1.index:
	this_seq=results1['sequence'][i]
	mobj=re.match("GT.*AG$",this_seq)
	if mobj:
		oseq=this_seq[0:10]+this_seq[len(this_seq)-10:len(this_seq)]
		if len(oseq)==20:
			print(oseq,file=out1)
	mobj2=re.match("CT.*AC$",this_seq)
	if mobj2:
		r_seq=Seq.reverse_complement(this_seq)
		oseq=r_seq[0:10]+r_seq[len(this_seq)-10:len(this_seq)]
		if len(oseq)==20:
			print(oseq,file=out1)


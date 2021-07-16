import MySQLdb
import re
import pandas
import rpy2
import math
import numpy
from Bio import SeqIO
from rpy2.robjects.vectors import IntVector, FloatVector,StrVector
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr, data
from rpy2.robjects import pandas2ri
pandas2ri.activate()
database_name='pse'
db=MySQLdb.connect("localhost", "root", "", database_name, charset='utf8' )
command_intron1="select * from chromosome_pflava where gene_number=1 and telomere=2 and intron_number=1;"
results1=pandas.read_sql(command_intron1,db)
command_intron2="select * from chromosome_pcarnea where gene_number=1 and telomere=2 and intron_number=1;"
results2=pandas.read_sql(command_intron2,db)
pos_list1=[0]*100
pos_list2=[0]*100
intron_count=0
for i in results1.index:
	num=results1['intron_number'][i]
	pos_str=results1['intron_pos_list'][i]
	f_length=results1['length'][i]
	f_strand=results1['single_gene_strand'][i]
	f_pos_list=pos_str.split(",")
	if f_strand=='+':
		for j in range(len(f_pos_list)):
			rel_pos=int(int(f_pos_list[j])/f_length*100)
			intron_count+=1
			pos_list1[rel_pos]+=1
	else:
		for j in range(len(f_pos_list)):
			rel_pos=100-int(int(f_pos_list[len(f_pos_list)-1-j])/f_length*100)
			intron_count+=1
			pos_list1[rel_pos]+=1
for i in range(len(pos_list1)):
	pos_list1[i]=pos_list1[i]/intron_count
yl1=ro.FloatVector(pos_list1)
xl1=ro.IntVector(numpy.arange(0,100))

intron_count=0
for i in results2.index:
	num=results2['intron_number'][i]
	pos_str=results2['intron_pos_list'][i]
	f_length=results2['length'][i]
	f_strand=results2['single_gene_strand'][i]
	f_pos_list=pos_str.split(",")
	if f_strand=='+':
		for j in range(len(f_pos_list)):
			rel_pos=int(int(f_pos_list[j])/f_length*100)
			intron_count+=1
			pos_list2[rel_pos]+=1
	else:
		for j in range(len(f_pos_list)):
			rel_pos=100-int(int(f_pos_list[len(f_pos_list)-1-j])/f_length*100)
			intron_count+=1
			pos_list2[rel_pos]+=1
for i in range(len(pos_list2)):
	pos_list2[i]=pos_list2[i]/intron_count
yl=ro.FloatVector(pos_list2)
xl=ro.IntVector(numpy.arange(0,100))
dfr=ro.r['data.frame'](x=xl,y=yl,z=yl1)
pp = ggplot2.ggplot(dfr) +\
	ggplot2.geom_point(ggplot2.aes_string(x='x',y='y'),color="red",size=3)+\
	ggplot2.geom_point(ggplot2.aes_string(x='x',y='z'),color="green",size=3)
pp.plot()
ggsave=ro.r["ggsave"]
ggsave(r"C:\Users\wbzhe\papers\pse\Figures\material\intron_distribution_pse.jpg",pp,width=12,height=4,dpi=600)
x=input()
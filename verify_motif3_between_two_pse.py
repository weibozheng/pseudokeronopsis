import MySQLdb
import re
import pandas
from sqlalchemy import create_engine
from Bio import SeqIO
import numpy
from math import isnan
import sys
database_name="pse"
table_name="chromosome_pcarnea"
db = MySQLdb.connect("localhost", "root", "", database_name, charset='utf8' )
sql_cmd="select * from "+table_name+";"
results=pandas.read_sql(sql_cmd,db)
count=0
dict_carnea=dict()
for i in results.index:
	tag=str(results['chromosome'][i])
	if not isnan(results['motif3_pos'][i]):
		dict_carnea.update({tag:1})
	else:
		dict_carnea.update({tag:0})
table_name="chromosome_pflava"
sql_cmd="select * from "+table_name+";"
results2=pandas.read_sql(sql_cmd,db)
count=0
dict_flava=dict()
for i in results2.index:
	tag=str(results2['chromosome'][i])
	if not isnan(results2['motif3_pos'][i]):
		dict_flava.update({tag:1})
	else:
		dict_flava.update({tag:0})
table_name="brh_chromosome"
sql_cmd="select * from "+table_name+";"
results3=pandas.read_sql(sql_cmd,db)

count_no=0
count_c=0
count_f=0
count_cf=0
count_pc=0
count_pf=0
for i in results3.index:
	tag1=results3['chromosome_pcarnea_chromosome'][i]
	tag2=results3['chromosome_pflava_chromosome'][i]
	if dict_carnea[tag1]==1 and dict_flava[tag2]==1:
		count_cf+=1
		count_pc+=1
	if dict_carnea[tag1]==0 and dict_flava[tag2]==0:
		count_no+=1
		count_pf+=1
	if dict_carnea[tag1]==1 and dict_flava[tag2]==0:
		count_c+=1
		count_pc+=1
		
	if dict_carnea[tag1]==0 and dict_flava[tag2]==1:
		count_f+=1
		count_pf+=1
print(count_no,count_c,count_f,count_cf,count_pc,count_pf)

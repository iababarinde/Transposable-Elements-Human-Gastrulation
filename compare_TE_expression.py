def compare_TE_expression(inmatrix,index1,index2,map_file,TE_file,coding_file,outfile,normalizer):
	import os
	type_dict=dict()
	subtype_dict=dict()
	type_map=dict()
	subtype_map=dict()
	te_map=dict()
	te_dict=dict()
	lc,pc=-1,0
	with open(map_file,'r') as yui:
		for line in yui:
			split=line.split()
			lc+=1
			if lc==0:
				continue
			te_tag=title_tag=split[3]+':'+split[4]+':'+split[0]
			subtype_tag=split[3]+':'+split[4]
			type_tag=split[3]
			subtype_map[split[0]]=subtype_tag
			subtype_dict[subtype_tag]=dict()
			type_map[split[0]]=type_tag
			type_dict[type_tag]=dict()
			te_map[split[0]]=te_tag
			te_dict[te_tag]=dict()
	coding_dict=dict()
	codmap={'0':'noncoding','1':'coding'}
	with open(coding_file,'r') as yui:
		for line in yui:
			split=line.split()
			if (split[-1] in codmap):
				coding_dict[split[0]]=codmap[split[-1]]
	print (len(coding_dict),'target transcripts found in the coding file')
	if '.gz' in TE_file:
		import gzip
		yui=gzip.open(TE_file,'r')
	else:
		yui=open(TE_file,'r')
	tid_dict=dict()
	for line in yui:
		if '.gz' in TE_file:
			line=line.decode()
		if line[0]=='#':
			continue
		split=line.split()
		if split[2] in te_map:
			tid=split[0].split('::')[0]
			te=split[2]
			if tid in coding_dict:
				codtye=coding_dict[tid]
				te_name,type_name,subtype_name=te_map[te],type_map[te],subtype_map[te]
				if te_name not in te_dict:
					te_dict[te_name]=dict()
				if type_name not in type_dict:
					type_dict[type_name]=dict()
				if subtype_name not in subtype_dict:
					subtype_dict[subtype_name]=dict()
				te_dict[te_name][tid],type_dict[type_name][tid],subtype_dict[subtype_name][tid]=0,0,0
				tid_dict[tid]=0
	yui.close()
	print (len(te_dict),'TEs found')
	headers=[]
	from numpy import log2,mean,median
	expr_dict=dict()
	raw_expr_dict=dict()
	lc=-1
	ind1,ind2=[],[]
	all_dict={'coding':[[],[]],'noncoding':[[],[]]}
	raw_all_dict={'coding':[[],[]],'noncoding':[[],[]]}
	with open(inmatrix,'r') as yui:
		for line in yui:
			lc+=1
			split=line.split()
			if lc==0:
				for one in set(index1):
					if one in split:
						ind1.append(split.index(one))
					elif isinstance(one,int) and one<len(split):
						ind1.append(one)
					else:
						try:
							ind1.append(int(one))
						except:
							return ('Something is wrong with index1:',one)
				for one in set(index2):
					if one in split:
						ind2.append(split.index(one))
					elif isinstance(one,int) and one<len(split):
						ind2.append(one)
					else:
						try:
							ind2.append(int(one))
						except:
							return ('Something is wrong with index2:',one)
			elif sum([float(i) for i in split[1:]])==0:
				continue
			elif (split[0] in tid_dict):
				val1=log2(mean([float(split[i]) for i in ind1])+normalizer)
				val2=log2(mean([float(split[i]) for i in ind2])+normalizer)
				raw1,raw2=mean([float(split[i]) for i in ind1]),mean([float(split[i]) for i in ind2])
				if raw1==raw2==0:
					continue
				all_dict[coding_dict[split[0]]][0].append(val1)
				all_dict[coding_dict[split[0]]][1].append(val2)
				expr_dict[split[0]]=[val1,val2]
				raw_all_dict[coding_dict[split[0]]][0].append(raw1)
				raw_all_dict[coding_dict[split[0]]][1].append(raw2)
				raw_expr_dict[split[0]]=[raw1,raw2]
	codcount,noncount=1,1
	for group in [type_dict,subtype_dict,te_dict]:
		for te in group:
			cdict={'coding':[[],[]],'noncoding':[[],[]]}
			for tid in group[te]:
				if tid in expr_dict:
					cod=coding_dict[tid]
					cdict[cod][0].append(expr_dict[tid][0])
					cdict[cod][1].append(expr_dict[tid][1])
			if len(cdict['coding'][0])>2:
				codcount+=1
			if len(cdict['noncoding'][0])>2:
				noncount+=1
	print ('Number of coding TE tests:',codcount)
	print ('Number of noncoding TE tests:',noncount)
	from scipy.stats import wilcoxon, ttest_rel
	newcod=open(outfile+'_coding','w')
	newnon=open(outfile+'_noncoding','w')
	newcod.write('TE\tcodtype\tnumber\traw_mean1\traw_mean2\tmean1\tmean2\tmedian1\tmedian2\tratio\tlfc\tttest_pvalue\tttest_adj_value\twilcoxon_pvalue\twilcoxon_adj_pvalue\n')
	newnon.write('TE\tcodtype\tnumber\traw_mean1\traw_mean2\tmean1\tmean2\tmedian1\tmedian2\tratio\tlfc\tttest_pvalue\tttest_adj_value\twilcoxon_pvalue\twilcoxon_adj_pvalue\n')
	ts,tp=ttest_rel(all_dict['coding'][0],all_dict['coding'][1])
	ws,wp=wilcoxon(all_dict['coding'][0],all_dict['coding'][1])
	utp="{:.2e}".format(tp)
	utq="{:.2e}".format(tp*(codcount))
	uwp="{:.2e}".format(wp)
	uwq="{:.2e}".format(wp*(codcount))
	mean1,mean2=mean(all_dict['coding'][0]),mean(all_dict['coding'][1])
	med1,med2=median(all_dict['coding'][0]),median(all_dict['coding'][1])
	rat=round(log2((mean1+normalizer)/(mean2+normalizer)),3)
	raw_mean1,raw_mean2=mean(raw_all_dict['coding'][0]),mean(raw_all_dict['coding'][1])
	lfc=round(log2((raw_mean1+normalizer)/(raw_mean2+normalizer)),3)
	newcod.write('allTE\tcoding\t'+str(len(all_dict['coding'][0]))+'	'+str(raw_mean1)+'	'+str(raw_mean2)+'	'+str(mean1)+'	'+str(mean2)+'	'+str(med1)+'	'+str(med2)+'	'+str(rat)+'	'+str(lfc)+'	'+utp+'	'+utq+'	'+uwp+'	'+uwq+'\n')
	tt,tp=ttest_rel(all_dict['noncoding'][0],all_dict['noncoding'][1])
	wt,wp=wilcoxon(all_dict['noncoding'][0],all_dict['noncoding'][1])
	utp="{:.2e}".format(tp)
	utq="{:.2e}".format(tp*(noncount))
	uwp="{:.2e}".format(wp)
	uwq="{:.2e}".format(wp*(noncount))
	mean1,mean2=mean(all_dict['noncoding'][0]),mean(all_dict['noncoding'][1])
	med1,med2=median(all_dict['noncoding'][0]),median(all_dict['noncoding'][1])
	rat=round(log2((mean1+normalizer)/(mean2+normalizer)),3)
	raw_mean1,raw_mean2=mean(raw_all_dict['noncoding'][0]),mean(raw_all_dict['noncoding'][1])
	lfc=round(log2((raw_mean1+normalizer)/(raw_mean2+normalizer)),3)
	newnon.write('allTE\tnoncoding\t'+str(len(all_dict['noncoding'][0]))+'	'+str(raw_mean1)+'	'+str(raw_mean2)+'	'+str(mean1)+'	'+str(mean2)+'	'+str(med1)+'	'+str(med2)+'	'+str(rat)+'	'+str(lfc)+'	'+utp+'	'+utq+'	'+uwp+'	'+uwq+'\n')
	dc,dn=0,0
	for group in [type_dict,subtype_dict,te_dict]:
		for te in group:
			cdict={'coding':[[],[]],'noncoding':[[],[]]}
			raw_cdict={'coding':[[],[]],'noncoding':[[],[]]}
			for tid in group[te]:
				if tid in expr_dict:
					cod=coding_dict[tid]
					cdict[cod][0].append(expr_dict[tid][0])
					cdict[cod][1].append(expr_dict[tid][1])
					raw_cdict[cod][0].append(raw_expr_dict[tid][0])
					raw_cdict[cod][1].append(raw_expr_dict[tid][1])
			if len(cdict['coding'][0])>2:
				dc+=1
				ts,tp=ttest_rel(cdict['coding'][0],cdict['coding'][1])
				try:
					ws,wp=wilcoxon(cdict['coding'][0],cdict['coding'][1])
				except:
					wp=1
				utp="{:.2e}".format(tp)
				utq="{:.2e}".format(min(1,tp*(codcount)))
				uwp="{:.2e}".format(wp)
				uwq="{:.2e}".format(min(1,wp*(codcount)))
				mean1,mean2=mean(cdict['coding'][0]),mean(cdict['coding'][1])
				med1,med2=median(cdict['coding'][0]),median(cdict['coding'][1])
				rat=round(log2((mean1+normalizer)/(mean2+normalizer)),3)
				raw_mean1,raw_mean2=mean(raw_cdict['coding'][0]),mean(raw_cdict['coding'][1])
				lfc=round(log2((raw_mean1+normalizer)/(raw_mean2+normalizer)),3)
				newcod.write(te+'\tcoding\t'+str(len(cdict['coding'][0]))+'	'+str(raw_mean1)+'	'+str(raw_mean2)+'	'+str(mean1)+'	'+str(mean2)+'	'+str(med1)+'	'+str(med2)+'	'+str(rat)+'	'+str(lfc)+'	'+utp+'	'+utq+'	'+uwp+'	'+uwq+'\n')
			if len(cdict['noncoding'][0])>2:
				dn+=1
				tt,tp=ttest_rel(cdict['noncoding'][0],cdict['noncoding'][1])
				try:
					wt,wp=wilcoxon(cdict['noncoding'][0],cdict['noncoding'][1])
				except:
					wp=1
				utp="{:.2e}".format(tp)
				utq="{:.2e}".format(min(1,tp*(noncount)))
				uwp="{:.2e}".format(wp)
				uwq="{:.2e}".format(min(1,wp*(noncount)))
				mean1,mean2=mean(cdict['noncoding'][0]),mean(cdict['noncoding'][1])
				med1,med2=median(cdict['noncoding'][0]),median(cdict['noncoding'][1])
				rat=round(log2((mean1+normalizer)/(mean2+normalizer)),3)
				raw_mean1,raw_mean2=mean(raw_cdict['noncoding'][0]),mean(raw_cdict['noncoding'][1])
				lfc=round(log2((raw_mean1+normalizer)/(raw_mean2+normalizer)),3)
				newnon.write(te+'\tnoncoding\t'+str(len(cdict['noncoding'][0]))+'	'+str(raw_mean1)+'	'+str(raw_mean2)+'	'+str(mean1)+'	'+str(mean2)+'	'+str(med1)+'	'+str(med2)+'	'+str(rat)+'	'+str(lfc)+'	'+utp+'	'+utq+'	'+uwp+'	'+uwq+'\n')
	print ('Coding TE computed:',dc)
	print ('Noncoding TE computed:',dn)
	newcod.close()
	newnon.close()


from optparse import OptionParser, OptionGroup
import os,sys,gzip
parser = OptionParser(usage="python3 %prog [options]", version="%prog version 1.0")

parser.add_option('-x','--expression', dest='expression', help='Absolute path to the expression file [required]',type='str')
parser.add_option('-c','--coding_file', dest='coding_file', help='Absolute path to the feelnc RF file [required]',type='str')
parser.add_option('-t','--te_file', dest='te_file', help='Absolute path to the TE file mapping TE to transcript [required]',type='str')
parser.add_option('-m','--map_file', dest='map_file', help='Absolute path to the TE annotation file, mapping TEs to type and subtype [required]',type='str')
parser.add_option('-o','--outfile', dest='outfile', help='Absolute path to file to write the final output [required]',type='str')
parser.add_option('-1','--index1', action='append', dest='index1', help='Column or column name for condition 1. Can be supplied multiple times [required]',type='str')
parser.add_option('-2','--index2', action='append', dest='index2', help='Column or column name for condition 2. Can be supplied multiple times [required]',type='str')
parser.add_option('-n','--normalizer', dest='normalizer', help='Value to add before log2 transformation [Default:1]',type='float', default=1)



(options,args)=parser.parse_args()


if (options.expression==None) or (options.te_file==None) or (options.map_file==None) or (options.coding_file==None) or (options.outfile==None):
	print ('\nRequired filed(s) not supplied\n#  This compares TE expressions\n')
	parser.print_help()
	sys.exit(1)

if (options.index1==None) or (options.index1==None):
	print ('\nRequired filed(s) not supplied\n# idex1 and index2 are required. This compares TE expressions\n')
	parser.print_help()
	sys.exit(1)


myarg=vars(options)

for arg in myarg:
	print (arg.upper(),':',myarg[arg])

compare_TE_expression(options.expression,options.index1,options.index2,options.map_file,options.te_file,options.coding_file,options.outfile,options.normalizer)
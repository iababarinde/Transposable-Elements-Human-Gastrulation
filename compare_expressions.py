def compare_expressions(filesA, filesB, coding_file, outfile, columnA, columnB, normalizer, xlabel, ylabel, main_title):
	if len(filesA)!=len(columnA) or len(filesB)!=len(columnB):
		print ('Files and columns inconsistencies, please check')
		return
	import os
	feelnc_dict=dict()
	codict={'1':'coding','0':'noncoding'}
	codcount={'coding':0,'noncoding':0}
	with open(coding_file,'r') as yui:
		for line in yui:
			split=line.split()
			try:
				feelnc_dict[split[0]]=codict[split[-1]]
				codcount[codict[split[-1]]]+=1
			except:
				pass
	print ('Feelnc Results:',codcount)
	base_dict=dict()
	print ('File(s) A...')
	fc=0
	union_dict=dict()
	for a in range(len(filesA)):
		fc+=1
		tc,nz=0,0
		one=filesA[a]
		with open(one,'r') as current:
			for line in current:
				split=line.split()
				try:
					value=float(split[columnA[a]])
					tc+=1
					if value>0:
						nz+=1
				except:
					if columnA[a] in split:
						columnA[a]=split.index(columnA[a])
					else:
						columnA[a]=int(columnA[a])
					continue
				try:
					base_dict[split[0].split('::')[0]].append(value)
				except:
					base_dict[split[0].split('::')[0]]=[value]
				union_dict[split[0].split('::')[0]]=0
		print (os.path.split(one)[-1],tc,nz)
	print ('Number of base files:',fc)
	print ('Union number of transcripts:',len(union_dict))
	focus_dict=dict()
	print ('File(s) B...')
	fc=0
	for b in range(len(filesB)):
		fc+=1
		tc,nz=0,0
		one=filesB[b]
		with open(one,'r') as current:
			for line in current:
				split=line.split()
				try:
					value=float(split[columnB[b]])
					tc+=1
					if value>0:
						nz+=1
				except:
					if columnB[b] in split:
						columnB[b]=split.index(columnB[b])
					else:
						columnB[b]=int(columnB[b])
					continue
				try:
					focus_dict[split[0].split('::')[0]].append(value)
				except:
					focus_dict[split[0].split('::')[0]]=[value]
				union_dict[split[0].split('::')[0]]=0
		print (os.path.split(one)[-1],tc,nz)
	print ('Number of base files:',fc)
	print ('Union number of transcripts:',len(union_dict))
	from numpy import mean, median,log2
	from scipy import stats
	cx,cy,nx,ny=[],[],[],[]
	uset=0
	for trans in union_dict:
		bvalue=0 if (trans not in base_dict) else mean(base_dict[trans])
		fvalue=0 if (trans not in focus_dict) else mean(focus_dict[trans])
		if bvalue==0 and fvalue==0:
			continue
		uset+=1
		coding='coding' if (trans in feelnc_dict and feelnc_dict[trans]=='coding') else 'noncoding'
		if coding=='coding':
			cx.append(log2(bvalue+normalizer))
			cy.append(log2(fvalue+normalizer))
		else:
			nx.append(log2(bvalue+normalizer))
			ny.append(log2(fvalue+normalizer))
	
	print ('Useful transcripts:',uset)
	print ('	coding:',len(cx))
	print ('	Noncoding:',len(nx))

	ccor,p=stats.spearmanr(cy,cx)
	ccor1,p1=stats.pearsonr(cy,cx)
	print ('Coding transcript correlations:')
	print ('	spearman:',ccor)
	print ('	pearson:',ccor1)


	ncor,p=stats.spearmanr(ny,nx)
	ncor1,p1=stats.pearsonr(ny,nx)
	print ('Noncoding transcript correlations:')
	print ('	spearman:',ncor)
	print ('	pearson:',ncor1)


	fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=False, figsize=(6, 2.5))

	hb = ax0.hexbin(cx, cy, gridsize=50, bins='log', cmap='inferno', mincnt=1)
	xlim=(int(round(min(cx)-0.5)),int(round(max(cx)+0.5)))
	ylim=(int(round(min(cy)-0.5)),int(round(max(cy)+0.5)))
	ax0.set(xlim=xlim, ylim=ylim)
	ax0.set_title(main_title+": Coding: " + str(round(ccor,3)), fontsize=6)
	ax0.set_ylabel(ylabel)
	ax0.set_xlabel(xlabel)
	cb = fig.colorbar(hb, ax=ax0, label='Counts')
	plt.tight_layout()

	hb = ax1.hexbin(nx, ny, gridsize=50, bins='log', cmap='inferno', mincnt=1)
	xlim=(int(round(min(nx)-0.5)),int(round(max(nx)+0.5)))
	ylim=(int(round(min(ny)-0.5)),int(round(max(ny)+0.5)))
	ax1.set(xlim=xlim, ylim=ylim)
	ax1.set_title(main_title+": Noncoding: " + str(round(ncor,3)), fontsize=6)
	cb = fig.colorbar(hb, ax=ax1, label='Counts')
	plt.tight_layout()
	fig.savefig(options.outfile)



import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from scipy import stats
import os,sys


parser = OptionParser()
parser.add_option('-a','--fileA', action='append', dest='fileA', help='Absolute path to the base (reference) gtf file [required]',type='str')
parser.add_option('-b','--fileB', action='append', dest='fileB', help='Absolute path to the focus (target) gtf file [required]',type='str')
parser.add_option('-c','--coding_file', dest='coding_file', help='Absolute path to the feelnc RF file [required]',type='str')
parser.add_option('-o','--outfile', dest='outfile', help='Absolute path to file to write the final output [required]',type='str')
parser.add_option('-1','--columnA', action='append', dest='columnA', help='Columns to extract in filesA, same number as filesA',type='str')
parser.add_option('-2','--columnB', action='append', dest='columnB', help='Columns to extract in filesB, same number as filesB',type='str')
parser.add_option('-n','--normalizer', dest='normalizer', help='Normalizer to add and to decide expression',type='float',default=0.001)
parser.add_option('-x','--xlab', dest='xlab', help='Labels for x axis of the plot',type='str', default='X-axis')
parser.add_option('-y','--ylab', dest='ylab', help='Labels for y axis of the plot',type='str', default='Y-axis')
parser.add_option('-m','--main', dest='main', help='Main title for the plot',type='str', default='')


(options,args)=parser.parse_args()


if (options.fileA==None) or (options.fileB==None) or (options.coding_file==None) or (options.outfile==None):
	print ('\nRequired filed(s) not supplied\n')
	parser.print_help()
	sys.exit(1)

myarg=vars(options)

for arg in myarg:
	print (arg.upper(),':',myarg[arg])


compare_expressions(options.fileA,options.fileB,options.coding_file,options.outfile,options.columnA,options.columnB,options.normalizer,options.xlab,options.ylab,options.main)
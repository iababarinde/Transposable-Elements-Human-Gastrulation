def splicing_TE_cov(gtf_list, TE_file_list, outfile, window, flank_count, map_file, coding_type='both', feelnc_file_list='',te_name='',te_subtype='',te_type='',color_list=None):
	import os
	value_dict=dict()
	map_dict=dict()
	lc,pc=-1,0
	with open(map_file,'r') as yui:
		for line in yui:
			split=line.split()
			lc+=1
			if lc==0:
				continue
			if te_name in split[:5]:
				map_dict[split[0]]=te_type+':'+te_subtype+':'+te_name
				title_tag=split[3]+':'+split[4]+':'+split[0]
				pc+=1
			elif te_subtype in split[:5]:
				map_dict[split[0]]=te_type+':'+te_subtype+':'+te_name
				title_tag=split[3]+':'+split[4]
				pc+=1
			elif te_type in split[:5]:
				map_dict[split[0]]=te_type+':'+te_subtype+':'+te_name
				title_tag=split[3]
				pc+=1
			elif te_name=='' and te_subtype=='' and te_type=='':
				map_dict[split[0]]=te_type+':'+te_subtype+':'+te_name
				title_tag='allTEs'
				pc+=1
	if pc==0:
		print ('Specified TE not found, exiting....')
		return
	coding_dict=dict()
	codmap={'1':'coding','0':'noncoding'}
	input_name_list=[]
	hold_value_dict=dict()
	color_dict=dict()
	for file_index in range(len(gtf_list)):
		gtf=gtf_list[file_index]
		TE_file=TE_file_list[file_index]
		name=os.path.split(gtf_list[file_index])[-1].replace('.gtf','').replace('.bed','')
		input_name_list.append(name)
		color_dict[name]=color_list[file_index]
		print (name,'...')
		if coding_type!='both':
			feelnc_file=feelnc_file_list[file_index]
			with open(feelnc_file,'r') as yui:
				for line in yui:
					split=line.split()
					if (split[-1] in codmap) and (codmap[split[-1]]==coding_type):
						coding_dict[split[0]]=0
			print ('\t',len(coding_dict),'transcripts found')
		my_dict=dict()
		if '.gtf' in gtf:
			with open(gtf,'r') as yui:
				for line in yui:
					split=line.split()
					if split[2]=='transcript':
						for i in range(6,len(split),1):
							if split[i]=='transcript_id':
								tid=split[i+1][1:-2]
								break
						my_dict[tid]=[split[6],[]]
					elif split[2]=='exon':
						coord=[int(split[3]),int(split[4])]
						coord.sort()
						my_dict[tid][-1].append(coord)
		elif '.bed' in gtf:
			with open(gtf,'r') as yui:
				for line in yui:
					split=line.split()
					tid=split[3]
					my_dict[tid]=[split[5],[[int(split[-1].split(',')[i])+int(split[1])+1,int(split[-1].split(',')[i])+int(split[-2].split(',')[i])+int(split[1])] for i in range(len(split[-1].split(',')))]]
		print ('\tTotal number of transcripts:',len(my_dict))
		junction_cov={0:0}
		junction_coord={0:[-window//2,window//2]}
		junction_list=[0]
		junction_tot={0:0}
		for i in range(1,flank_count+1,1):
			junction_list.extend([-i,i])
			junction_cov[-i]=0
			junction_cov[i]=0
			junction_tot[-i]=0
			junction_tot[i]=0
			junction_coord[-i]=[-i*window,(1-i)*window]
			junction_coord[i]=[(i-1)*window,i*window]
		junctions=dict()
		transcript_lengths=dict()
		junccount=0
		for tid in my_dict:
			if coding_type!='both' and (tid not in coding_dict):
				continue
			exons=my_dict[tid][-1]
			exons.sort()
			if len(exons)==1:
				continue
			exon_cov=[abs(one[1]-one[0])+1 for one in exons]
			junccount+=len(exon_cov)
			transcript_lengths[tid]=sum(exon_cov)
			if my_dict[tid][0]=='+':
				junctions[tid]=[sum(exon_cov[:i+1]) for i in range(len(exon_cov)-1)]
			else:
				junctions[tid]=[sum(exon_cov[::-1][:i+1]) for i in range(len(exon_cov)-1)]
			junction_tot[0]+=len(exon_cov)
			for juncs in junctions[tid]:
				for i in range(1,flank_count+1,1):
					if ((-i)*window)+juncs>=0:
						junction_tot[-i]+=1
					if ((i+1)*window)+juncs<=sum(exon_cov):
						junction_tot[i]+=1
		print ('\tTranscripts with junctions (multi-exon):',len(junctions))
		print ('\tTotal junction count:',junccount)
		print ('\tJunction bins count:',junction_tot)
		if '.gz' in TE_file:
			import gzip
			yui=gzip.open(TE_file,'r')
		else:
			yui=open(TE_file,'r')
		tid_te=dict()
		for line in yui:
			if '.gz' in TE_file:
				line=line.decode()
			if line[0]=='#':
				continue
			split=line.split()
			if split[2] in map_dict:
				tid=split[0].split('::')[0]
				if tid not in junctions:
					continue
				if coding_type!='both' and (tid not in coding_dict):
					continue
				acoord=[int(split[6]),int(split[7])]
				acoord.sort()
				if split[0] not in tid_te:
					tid_te[tid]=[]
				tid_te[tid].append(acoord)
		for tid in tid_te:
			pos=1
			current=tid_te[tid]
			current.sort()
			while pos<len(current):
				if current[pos-1][1]>=current[pos][0]:
					current[pos-1][1]=max(current[pos-1][1],current[pos][1])
					del current[pos-1]
				else:
					pos+=1
			for acoord in current:
				for junc in junctions[tid]:
					for j in junction_coord:
						if junction_coord[j][0]+junc<1 or junction_coord[j][1]+junc>transcript_lengths[tid]:
							continue
						if junction_coord[j][1]+junc>=acoord[0] and junction_coord[j][0]<=acoord[1]+junc:
							junction_cov[j]+=1		
		yui.close()
		junction_list.sort()
		hold_value=[]
		if file_index==0:
			for j in junction_list:
				value_dict[j]=dict()
		for j in junction_list:
			value_dict[j][name]=float(junction_cov[j])/junction_tot[j]
			hold_value.append(float(junction_cov[j])/junction_tot[j])
		hold_value_dict[name]=hold_value
	import matplotlib.pyplot as plt
	import numpy as np
	axis1 = plt.subplot(221)
	for name in input_name_list:
		axis1.plot(range(len(hold_value_dict[name])), hold_value_dict[name], label = name, linewidth = 1, color = color_dict[name])
	axis1.legend(bbox_to_anchor=(0.5, 1.2), loc='best',ncol=1)
	axis1.set_title(title_tag)
	#plt.tick_params('x', labelbottom=False, bottom=False)
	plt.xticks([0,float(len(junction_list)/2)-0.5,len(junction_list)-1],[junction_list[0]*window,0,junction_list[-1]*window],size=7)
	plt.axvline(x = float(len(junction_list)/2)-0.5, color = 'black', linestyle = '--', linewidth = 1)
	plt.tight_layout()
	plt.savefig(outfile+'.pdf')

	with open(outfile+'.tsv','w') as new:
		new.write('bin')
		for name in input_name_list:
			new.write('	'+name)
		new.write('\n')
		for j in junction_list:
			new.write(str(j))
			for name in input_name_list:
				new.write('	'+str(value_dict[j][name]))
			new.write('\n')
			


from optparse import OptionParser, OptionGroup
import os,sys,gzip
parser = OptionParser(usage="python3 %prog [options]", version="%prog version 1.0")

parser.add_option('-g', '--gtf', dest='gtf', action='append', help='Absolute path to the transcript file in gtf or bed format [required]',type='str')
parser.add_option('-c', '--coding_file', dest='coding_file', action='append', help='Absolute path to the feelnc RF file [required]',type='str')
parser.add_option('-t', '--te_file', dest='te_file', action='append', help='Absolute path to the TE file mapping TE to transcript [required]',type='str')
parser.add_option('-m', '--map_file', dest='map_file', help='Absolute path to the TE annotation file, mapping TEs to type and subtype [required]',type='str')
parser.add_option('-o', '--outfile', dest='outfile', help='Absolute path to file to write the final output [required]',type='str')
parser.add_option('--flank_count', dest='flank_count', help='Number of bins to extract from each side [default:10]',type='int',default=10)
parser.add_option('-w', '--window_size', dest='window_size', help='The size of a bin [default:10]',type='int',default=10)
parser.add_option('--te_name', dest='te_name', help='TE name to extract, e.g. AluYa5 [default:'' means everything]',type='str',default='')
parser.add_option('--te_subtype', dest='te_subtype', help='TE subtype to extract, e.g. L1 [default:'' means everything]',type='str',default='')
parser.add_option('--te_type', dest='te_type', help='TE type to extract, e.g. LINE [default:'' means everything]',type='str',default='')
parser.add_option('--codtype', dest='codtype', help='The coding type of the transcripts to consider[Options : coding, noncoding, both (default)]',type='str',default='both')
parser.add_option('--color_list', dest='color_list', action='append', help='color list for each sample [required] -i,-c,t and -n must be supplied in the same order',type='str')

(options,args)=parser.parse_args()


if (options.gtf==None) or (options.te_file==None) or (options.map_file==None) or (options.coding_file==None) or (options.outfile==None):
	print ('\nRequired filed(s) not supplied\n# This makes TE pileup for transcript regions\n')
	parser.print_help()
	sys.exit(1)

if not len(options.te_file)==len(options.coding_file):
	print ('\nNumbers of input gtf file, coding file and te file must be the same. This is not fulfiled\n')
	parser.print_help()
	sys.exit(2)

if (options.color_list==None):
	options.color_list=['blue','orange','green','red','purple','brown','pink','gray','yellow','cyan'][:len(options.te_file)]

myarg=vars(options)

for arg in myarg:
	print (arg.upper(),':',myarg[arg])


splicing_TE_cov(options.gtf, options.te_file, options.outfile, options.window_size, options.flank_count, options.map_file, options.codtype, options.coding_file,options.te_name,options.te_subtype,options.te_type,options.color_list)
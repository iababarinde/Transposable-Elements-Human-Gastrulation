def make_bins(start,end,number):
	pos=start
	my_list=[]
	for i in range(number):
		add=1 if end%number>i else 0
		npos=pos+(end//number)+add
		#print (pos,min(npos,end)-1,range(pos,npos), len(range(pos,npos)))
		my_list.append([pos,min(npos,end)-1])
		pos=npos
	return (my_list)


def noncoding_TE_overlap(gtf_file,feelnc,te_file,names_list,map_file,outfile,bcount,te_name,te_subtype,te_type, stranded_transcripts, color_list):
	import os
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
	print (lc,'lines found in the map file')
	print (pc,'TEs to be focussed')
	all_bin_dict=dict()
	color_dict=dict()
	for one_index in range(len(gtf_file)):
		if names_list==None:
			cname=os.path.split(gtf_file[one_index])[-1]
		else:
			cname=names_list[one_index]
		color_dict[cname]=color_list[one_index]
		print ('\n\nChecking',cname,'*********************************************')
		tc,uc=0,0
		te_dict=dict()
		if '.gz' in te_file[one_index]:
			yui=gzip.open(te_file[one_index],'r')
		else:
			yui=open(te_file[one_index],'r')
		for line in yui:
			line=line.decode()
			if line[0]=='#':
				continue
			tc+=1
			split=line.split()
			if split[2] in map_dict:
				chro=split[0].split(':')[0]
				acoord=[int(split[6]),int(split[7])]
				acoord.sort()
				tid=split[0].split('::')[0]
				if tid not in te_dict:
					te_dict[tid]=[]
				te_dict[tid].append(acoord)
				uc+=1
		yui.close()
		print (tc,'TE alignments in the TE file')
		print (uc,'useful TE alignments')
		codict=dict()
		cmap={'0':'noncoding','1':'coding'}
		count_dict={'1':0,'0':0}
		yui=open(feelnc[one_index],'r')
		for line in yui:
			split=line.split()
			try:
				codict[split[0]]=cmap[split[-1]]
				count_dict[split[-1]]+=1
			except:
				continue
		print (count_dict)
		bdict=dict()
		for one in range(bcount):
			bdict[one]=0
		tc,uc=0,0
		tbc,tebc=0,0
		with open(gtf_file[one_index],'r') as yui:
			tidict=dict()
			for line in yui:
				split=line.split()
				if len(split)<6 or line[0]=='#':
					continue
				if split[2]=='transcript':
					tc+=1
					for i in range(6,len(split),1):
						if split[i]=='transcript_id':
							tid=split[1+i][1:-2]
							tidict[tid]=[split[6],0]
							break
				elif split[2]=='exon':
					tidict[tid][1]+=abs(int(split[4])-int(split[3]))+1
			print ('Last ID:', tidict[tid])
			for tid in tidict:
				strand=tidict[tid][0]
				if (tid in codict) and (codict[tid]=='noncoding'):
					if tidict[tid][1]>=bcount:
						uc+=1
						c_bins=make_bins(0,tidict[tid][1],bcount)
						this=[0 for i in range(bcount)]
						if tid in te_dict:
							tbc+=1
							for i in range(bcount):
								cur=c_bins[i]
								for tcoord in te_dict[tid]:
									if tcoord[1]<cur[0]:
										continue
									elif tcoord[0]>cur[1]:
										break
									this[i]=1
									tebc+=1
									break
						if tidict[tid][0]=='+' or stranded_transcripts[0].lower() in 'ty':
							for i in range(bcount):
								if this[i]==1:
									bdict[i]+=1
						elif tidict[tid][0]=='-':
							this=this[::-1]
							for i in range(bcount):
								if this[i]==1:
									bdict[i]+=1
		print ('Total number of transdecoder predictions:',tc)
		print ('Number of useful transcripts:',uc)
		print ('Number of TE-containing  transcripts:',tbc)
		print ('Number of TE-containing bins:',tebc)
		all_bin_dict[cname]=[]
		for i in range(bcount):
			all_bin_dict[cname].append(float(bdict[i])/uc)
	import matplotlib.pyplot as plt
	import numpy as np
	fig, axs = plt.subplots(1, 1, figsize=(1.0, 1.4))
	#print (axs)
	axis1=axs
	for name in all_bin_dict:
		axis1.plot(range(len(all_bin_dict[cname])), all_bin_dict[name], label = name, linewidth=1, color = color_dict[name])
	#axis1.legend(bbox_to_anchor=(0.5, 1.2), loc='best',ncol=1)
	axis1.set_title(title_tag,fontsize = 8, pad='1.0')
	#plt.tick_params('x', labelbottom=False, bottom=False)
	plt.xticks([0,len(all_bin_dict[cname])-1], ['TSS','TES'],size=7)
	plt.yticks(size=7)
	axis1.tick_params(axis='x', pad=1)
	axis1.tick_params(axis='y', pad=1)
	plt.tight_layout()
	plt.savefig(outfile+'.pdf')
	name_list=list(all_bin_dict)
	with open(outfile+'.tsv','w') as new:
		new.write('bin')
		for name in name_list:
			new.write('	'+name)
		new.write('\n')
		for i in range(bcount):
			new.write(str(i))
			for name in name_list:
				new.write('	'+str(all_bin_dict[name][i]))
			new.write('\n')
 

	


from optparse import OptionParser, OptionGroup
import os,sys,gzip
parser = OptionParser(usage="python3 %prog [options]", version="%prog version 1.0")

parser.add_option('-i','--infile', dest='infile', action='append',  help='Absolute path to the assembly gtf file [required]',type='str')
parser.add_option('-c','--coding_file', dest='coding_file', action='append', help='Absolute path to the feelnc RF file [required]',type='str')
parser.add_option('-t','--te_file', dest='te_file', action='append', help='Absolute path to the TE file mapping TE to transcript [required]',type='str')
parser.add_option('-n','--names', dest='names', action='append', help='name for each sample [required] -i,-c,t and -n must be supplied in the same order',type='str')
parser.add_option('-m','--map_file', dest='map_file', help='Absolute path to the TE annotation file, mapping TEs to type and subtype [required]',type='str')
parser.add_option('-o','--outfile', dest='outfile', help='Absolute path to file to write the final output [required]',type='str')
parser.add_option('--bin_count', dest='bin_count', help='Number of bins to extract [default:10]',type='int',default=10)
parser.add_option('--te_name', dest='te_name', help='TE name to extract, e.g. AluYa5 [default:'' means everything]',type='str',default='')
parser.add_option('--te_subtype', dest='te_subtype', help='TE subtype to extract, e.g. L1 [default:'' means everything]',type='str',default='')
parser.add_option('--te_type', dest='te_type', help='TE type to extract, e.g. LINE [default:'' means everything]',type='str',default='')
parser.add_option('--stranded_transcripts', dest='stranded_transcripts', help='If the sequences of the stranded transcripts were used for the hmmer search [Options: yes, no, true, false]',type='str',default='no')
parser.add_option('--color_list', dest='color_list', action='append', help='color list for each sample [required] -i,-c,t and -n must be supplied in the same order',type='str')


(options,args)=parser.parse_args()


if (options.infile==None) or (options.te_file==None) or (options.map_file==None) or (options.coding_file==None) or (options.outfile==None):
	print ('\nRequired filed(s) not supplied\n# This makes TE pileup for noncoding transcripts\n')
	parser.print_help()
	sys.exit(1)

if not len(options.infile)==len(options.te_file)==len(options.coding_file):
	print ('\nNumbers of input gtf file, coding file and te file must be the same. This is not fulfiled\n')
	parser.print_help()
	sys.exit(2)


myarg=vars(options)

for arg in myarg:
	print (arg.upper(),':',myarg[arg])

if (options.color_list==None):
	options.color_list=['blue','orange','green','red','purple','brown','pink','gray','yellow','cyan'][:len(options.te_file)]

noncoding_TE_overlap(options.infile,options.coding_file,options.te_file,options.names,options.map_file,options.outfile,options.bin_count,options.te_name,options.te_subtype,options.te_type,options.stranded_transcripts,options.color_list)
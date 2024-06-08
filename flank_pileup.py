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


def region_TE_overlap(gtf,feelnc,te_file,map_file,outfile,bcount,te_name,te_subtype,te_type, codtype, region_length,side, center, strand, color_list, names_list):
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
	strand_map={'plus':'+','minus':'-'}
	strand=strand.lower()
	c_bins=make_bins(0,region_length,bcount)
	from numpy import mean
	color_dict=dict()
	for one_index in range(len(feelnc)):
		if len(feelnc)!=names_list:
			cname=os.path.split(gtf[one_index])[-1].replace('.txt','')
		else:
			cname=names_list[one_index]
		color_dict[cname]=color_list[one_index]
		print ('\n\nChecking',cname,'*********************************************')
		gtf_dict=dict()
		ttc,utc=0,0
		with open(gtf[one_index],'r') as yui:
			for line in yui:
				split=line.split()
				if len(split)>7 and split[2]=='transcript':
					ttc+=1
					if strand!='both' and strand_map[strand]!=split[6]:
						continue
					utc+=1
					for i in range(6,len(split),1):
						if split[i]=='transcript_id':
							tid=split[i+1][1:-2]
							gtf_dict[tid]=abs(int(split[3])-int(split[4]))+1
		print (ttc,'transcripts in the gtf file')
		print (utc,'transcripts with the required strand')
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
			if split[0].split('::')[0] not in gtf_dict:
				continue
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
			if split[0] not in gtf_dict:
				continue
			try:
				codict[split[0]]=cmap[split[-1]]
				count_dict[split[-1]]+=1
			except:
				continue
		print (count_dict)
		bdict=dict()
		denum_count=[]
		for one in range(bcount):
			bdict[one]=0
			denum_count.append(0)
		tc,uc,tbc,tebc=0,0,0,0
		for tid in codict:
			tc+=1
			if (codict[tid]==codtype) or (codtype=='both'):
				uc+=1
				this=[0 for i in range(bcount)]
				if tid in te_dict:
					tbc+=1
					for i in range(bcount):
						cur=c_bins[i]
						if side[:2]=='up' and (gtf_dict[tid]+region_length//2)<mean(cur):
							break
						if side[:4]=='down' and (region_length//2 - gtf_dict[tid])>mean(cur):
							break
						denum_count[i]+=1
						for tcoord in te_dict[tid]:
							if tcoord[1]<cur[0]:
								continue
							elif tcoord[0]>cur[1]:
								break
							this[i]=1
							tebc+=1
							break
				for i in range(bcount):
					if this[i]==1:
						bdict[i]+=1
		print ('Total number of transdecoder predictions:',tc)
		print ('Number of useful transcripts:',uc)
		print ('Bin count range:',min(denum_count),max(denum_count))
		print ('Number of TE-containing  transcripts:',tbc)
		print ('Number of TE-containing bins:',tebc)
		all_bin_dict[cname]=[]
		for i in range(bcount):
			all_bin_dict[cname].append(float(bdict[i])/denum_count[i])
	import matplotlib.pyplot as plt
	import numpy as np
	fig, axis1 = plt.subplots(1, 1, figsize=(1.4, 1.1))
	for name in all_bin_dict:
		axis1.plot(range(len(all_bin_dict[cname])), all_bin_dict[name], label = name, linewidth = 1, color = color_dict[name])
	#axis1.legend(bbox_to_anchor=(0.5, 1.2), loc='best',ncol=1)
	axis1.set_title(title_tag,fontsize = 8, pad='1.0')
	#plt.tick_params('x', labelbottom=False, bottom=False)
	plt.xticks([0,len(all_bin_dict[cname])/2.0,len(all_bin_dict[cname])-1],[-(region_length//2),center,region_length//2],size=7)
	plt.yticks(size=7)
	plt.axvline(x = len(all_bin_dict[cname])/2.0, color = 'black', linestyle = '--', linewidth = 1)
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

parser.add_option('-g','--gtf_file', dest='gtf_file', action='append', help='Absolute path to the gtf files [required]',type='str')
parser.add_option('-c','--coding_file', dest='coding_file', action='append', help='Absolute path to the feelnc RF file [required]',type='str')
parser.add_option('-t','--te_file', dest='te_file', action='append', help='Absolute path to the TE file mapping TE to transcript [required]',type='str')
parser.add_option('-m','--map_file', dest='map_file', help='Absolute path to the TE annotation file, mapping TEs to type and subtype [required]',type='str')
parser.add_option('-o','--outfile', dest='outfile', help='Absolute path to file to write the final output [required]',type='str')
parser.add_option('--bin_count', dest='bin_count', help='Number of bins to extract [default:10]',type='int',default=10)
parser.add_option('--te_name', dest='te_name', help='TE name to extract, e.g. AluYa5 [default:'' means everything]',type='str',default='')
parser.add_option('--te_subtype', dest='te_subtype', help='TE subtype to extract, e.g. L1 [default:'' means everything]',type='str',default='')
parser.add_option('--te_type', dest='te_type', help='TE type to extract, e.g. LINE [default:'' means everything]',type='str',default='')
parser.add_option('--codtype', dest='codtype', help='The coding type of the transcripts to consider[Options : coding, noncoding, both (default)]',type='str',default='both')
parser.add_option('--region_length', dest='region_length', help='The length of the region to consider [required]',type='int')
parser.add_option('--side', dest='side', help='The side of the transcripts being ckecked [options: up(stream) or down(stream)]',type='str')
parser.add_option('--center', dest='center', help='The label of the center [default: 0]',type='str',default='0')
parser.add_option('--strand', dest='strand', help='Strand of the transcript [Options: plus/minus/both, default:both]',type='str',default='both')
parser.add_option('--color_list', dest='color_list', action='append', help='color list for each sample [required] -g,-c,t and -n must be supplied in the same order',type='str')
parser.add_option('-n','--names', dest='names', action='append', help='name for each sample [required] -g,-c,t and -n must be supplied in the same order',type='str')



(options,args)=parser.parse_args()


if (options.te_file==None) or (options.map_file==None) or (options.coding_file==None) or (options.outfile==None):
	print ('\nRequired filed(s) not supplied\n# This makes TE pileup for transcript regions\n')
	parser.print_help()
	sys.exit(1)

if not len(options.gtf_file)==len(options.te_file)==len(options.coding_file):
	print ('\nNumbers of input gtf file, coding file and te file must be the same. This is not fulfiled\n')
	parser.print_help()
	sys.exit(2)

if (options.region_length==None):
	print ('--region_length is a required option. It MUST be supplied\n')
	parser.print_help()
	sys.exit(3)

if (options.color_list==None):
	options.color_list=['blue','orange','green','red','purple','brown','pink','gray','yellow','cyan'][:len(options.te_file)]

myarg=vars(options)

for arg in myarg:
	print (arg.upper(),':',myarg[arg])

region_TE_overlap(options.gtf_file,options.coding_file,options.te_file,options.map_file,options.outfile,options.bin_count,options.te_name,options.te_subtype,options.te_type,options.codtype,options.region_length,options.side,options.center,options.strand, options.color_list, options.names)

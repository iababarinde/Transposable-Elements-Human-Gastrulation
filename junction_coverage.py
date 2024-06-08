def junction_reads(gtf,alignments,outfile,anchor_length):
	import pysam, os
	my_dict=dict()
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
	print ('Total number of transcripts:',len(my_dict))
	junctions=dict()
	for tid in my_dict:
		exons=my_dict[tid][-1]
		exons.sort()
		if len(exons)==1:
			continue
		exon_cov=[abs(one[1]-one[0])+1 for one in exons]
		junclist=[]
		if my_dict[tid][0]=='+':
			junctions[tid]=[[sum(exon_cov[:i+1]),0] for i in range(len(exon_cov)-1)]
		else:
			junctions[tid]=[[sum(exon_cov[::-1][:i+1]),0] for i in range(len(exon_cov)-1)]
	print ('Transcripts with junctions (multi-exon):',len(junctions))
	ar,jr=0,0
	for alignment in alignments:
		ar0,jr0=0,0
		with pysam.Samfile(alignment, 'rb') as samfile:
			for read in samfile.fetch(until_eof=True): #samfile.fetch(until_eof=True):
				ar0+=1
				if read.is_unmapped: # or (read.opt('NH')>1 and (read.qname in all_dict)):
					continue
				coord=read.positions
				coord=[coord[0],coord[-1]]
				coord.sort()
				tid=read.reference_name.split('::')[0]
				if tid not in junctions:
					continue
				try:
					found=0
					for i in range(len(junctions[tid])):
						if coord[0]+anchor_length<=junctions[tid][i][0]<=coord[1]-anchor_length:
							junctions[tid][i][1]+=1
							jr0+=1
							found=1	
				except:
					pass
		print (os.path.split(alignment)[-1],ar0,jr0)
		ar+=ar0
		jr+=jr0
	print ('Total number of reads:',ar)
	print ('Total number of junction reads:',jr)
	tjc,cjc=0,0
	with open(outfile,'w') as new:
		for tid in junctions:
			for i in range(len(junctions[tid])):
				new.write(tid+'	'+str(i)+'	'+str(junctions[tid][i][1])+'\n')
				tjc+=1
				if junctions[tid][i][1]>0:
					cjc+=1
	print ('Total junction count:',tjc)
	print ('Covered junction count:',cjc)

#junction_reads('/data3/isaac/early_diff/NEW/filtered/H1.gtf','/data3/isaac/early_diff/NEW/quantification/bowtie2/H1_btw2_rp1.bam','/data3/isaac/early_diff/NEW/quantification/bowtie2/H1_btw2_rp1.junctions')

from optparse import OptionParser, OptionGroup
import os,sys
parser = OptionParser(usage="python3 %prog [options]", version="%prog version 1.0")

parser.add_option('-g','--gtf_file', dest='gtf_file', help='Absolute path to the gtf file [required]',type='str')
parser.add_option('-b','--bam_file', dest='bam_file',action='append', help='Absolute path to the alignment bam file [required], can be supplied multiple times',type='str')
parser.add_option('-o','--outfile', dest='outfile', help='Absolute path to file to write the final output [required]',type='str')
parser.add_option('-a','--anchor_length', dest='anchor_length', help='Anchor length to be satisfied at the two ends of the junction [default=0, i.e no filter]',type='int',default=0)


(options,args)=parser.parse_args()


if (options.gtf_file==None) or (options.bam_file==None) or (options.outfile==None):
	print ('\nRequired filed(s) not supplied\n# #This produces the \n')
	parser.print_help()
	sys.exit(1)

myarg=vars(options)

for arg in myarg:
	print (arg.upper(),':',myarg[arg])

junction_reads(options.gtf_file,options.bam_file,options.outfile,options.anchor_length)
def transcript_overlap(infile_list,outfile, min_cov=0, max_hit=1000):
	import os
	#print ('*****************',os.path.split(infile)[-1])
	rdict=dict()
	ta,cov_pass=0,0
	for infile in infile_list:
		cta,ccov_pass=0,0
		with open(infile,'r') as yui:
			for line in yui:
				cta+=1
				split=line.split()
				cov=float(abs(int(split[3])-int(split[2])))*100/int(split[1])
				if cov<min_cov:
					continue
				nmb=int(split[9])
				trans_nmb=float(nmb)/int(split[6])
				ccov_pass+=1
				if (split[0] not in rdict) or ([nmb,trans_nmb,cov]>rdict[split[0]][0][:3]):
					rdict[split[0]]=[[nmb,trans_nmb,cov,split[5]]]
				elif [nmb,trans_nmb,cov]==rdict[split[0]][0][:3]:
					rdict[split[0]].append([nmb,trans_nmb,cov,split[5]])
		print ('::::',os.path.split(infile)[-1],'::::',cta,ccov_pass)
		ta+=cta
		cov_pass+=ccov_pass
	print ('Total alignment:',ta)
	print ('Unique reads:',len(rdict))
	print ('Reads passing min_cov filter:',cov_pass)
	count_dict=dict()
	mhp=0
	unique,multi=0,0
	for read in rdict:
		if len(rdict[read])==1:
			unique+=1
		else:
			multi+=1
		if len(rdict[read])>max_hit:
			continue
		mhp+=1
		for tid in rdict[read]:
			if tid[-1] not in count_dict:
				count_dict[tid[-1]]=0
			count_dict[tid[-1]]+=1
	print ('Number of unique mapping reads:',unique)
	print ('Number of multi-mapping reads:',multi)
	print ('Reads passing max_hit filter:',mhp)
	print ('Number of transcripts:',len(count_dict))
	with open(outfile,'a') as new:
		for tid in count_dict:
			new.write(tid.split('::')[0]+'	'+str(count_dict[tid])+'	'+str(float(count_dict[tid])*1000000/mhp)+'\n')


from optparse import OptionParser, OptionGroup
import os,sys
parser = OptionParser(usage="python3 %prog [options]", version="%prog version 1.0")

parser.add_option('-i','--infile', dest='infile', action='append', help='Absolute path to the minimap paf alignment output, can be supplied multiple times [required]',type='str')
parser.add_option('-o','--outfile', dest='outfile', help='Absolute path to file to write the final output [required]',type='str')
parser.add_option('--min_cov', dest='min_cov', help='Minimum read coverage. Reads with less coverage are discarded [default:0]',type='int', default=0)
parser.add_option('--max_hit', dest='max_hit', help='Maximum number of hits to consider. Reads with more reads are discarded [default:1000]',type='int', default=1000)



(options,args)=parser.parse_args()


if (options.infile==None) or (options.outfile==None):
	print ('\nRequired filed(s) not supplied\n# #This produces the \n')
	parser.print_help()
	sys.exit(1)

myarg=vars(options)

for arg in myarg:
	print (arg.upper(),':',myarg[arg])

transcript_overlap(options.infile,options.outfile,options.min_cov,options.max_hit)
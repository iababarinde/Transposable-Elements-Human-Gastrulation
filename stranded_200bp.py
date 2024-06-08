def stranded_200bp(infile,outfile):
	itc,iec=0,0
	wtc,wec=0,0
	cec=0
	new=open(outfile,'w')
	tel=0
	hold=''
	with open(infile,'r') as yui:
		for line in yui:
			split=line.split()
			if line[0]=='#' or len(split)<7 or (split[6] not in '+-'):
				continue
			if split[2]=='transcript':
				if tel>=200:
					new.write(hold)
					wtc+=1
					wec+=cec
				cec,tel,hold=0,0,''
				itc+=1
			elif split[2]=='exon':
				iec+=1
				cec+=1
				tel+=abs(int(split[4])-int(split[3]))+1
			hold+=line
	if tel>=200:
		new.write(hold)
		wtc+=1
		wec+=cec
	new.close()
	import os
	print ('File:',os.path.split(infile)[-1])
	print ('\tInput transcript count:',itc)
	print ('\tInput exon count:',iec)
	print ('\tInput exon per transcript:',float(iec)/itc)
	print ('\tFiltered transcript count:',wtc)
	print ('\tFiltered exon count:',wec)
	print ('\tFiltered exon per transcript:',float(wec)/wtc)

import os
for one in os.listdir('/data3/isaac/early_diff/NEW/stringtie'):
	stranded_200bp(os.path.join('/data3/isaac/early_diff/NEW/stringtie',one),os.path.join('/data3/isaac/early_diff/NEW/filtered',one))
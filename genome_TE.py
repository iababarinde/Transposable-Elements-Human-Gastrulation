def extract_coverage(infile,outfile, genome_size):
	my_dict=dict()
	with open(infile,'r') as yui:
		for line in yui:
			split=line.split()
			te_type=split[3].split(':')[-1]
			coord=[int(split[1]),int(split[2])]
			coord.sort()
			if te_type not in my_dict:
				my_dict[te_type]=dict()
			if split[0] not in my_dict[te_type]:
				my_dict[te_type][split[0]]=[]
			my_dict[te_type][split[0]].append(coord)
	count_dict=dict()
	for te_type in my_dict:
		tot=0
		bc,ac=0,0
		for chro in my_dict[te_type]:
			current=my_dict[te_type][chro]
			current.sort()
			pos=1
			bc+=len(current)
			while pos<len(current)-1:
				if current[pos][0]<=current[pos-1][1]:
					current[pos-1][1]=max(current[pos-1][1],current[pos][1])
					del current[pos]
				else:
					tot+=current[pos-1][1]-current[pos-1][0]+1
					pos+=1
			tot+=current[pos-1][1]-current[pos-1][0]+1
			ac+=len(current)
		print (te_type, bc, ac, tot)
		count_dict[te_type]=tot
	perc_dict=dict()
	with open(outfile+'_TE_count.tsv','w') as new:
		for te in my_dict:
			print (te, count_dict[te], count_dict[te]*100.0/genome_size)
			perc_dict[te]=count_dict[te]*100.0/genome_size
			new.write(te+'	'+str(count_dict[te])+'	'+str(count_dict[te]*100.0/genome_size)+'\n')
	#focus_te=['SINE', 'Retroposon', 'LINE', 'DNA', 'LTR', 'snRNA']
	focus_te=['LINE', 'SINE', 'LTR', 'DNA', 'Retroposon', 'snRNA']
	perc_list=[perc_dict[te] for te in focus_te]
	count_list=[count_dict[te] for te in focus_te]
	round_perc=[round(perc_dict[te],2) for te in focus_te]
	from matplotlib import pyplot as plt
	import numpy as np
	import matplotlib as mpl
	ax = plt.subplot(361)
	p1=ax.bar(focus_te, perc_list, color='grey')
	ax.bar_label(p1, labels=round_perc, label_type='center', size=7, rotation=90)
	plt.xticks(size=7, rotation =90)
	ax.set_title('Genome', fontsize = 8)
	ax.set_ylabel('Coverage (%)', fontsize = 8)
	plt.yticks(size=7)
	plt.xticks(size=7)
	plt.tight_layout()
	plt.savefig(outfile+'_barplots.pdf')
	plt.close()


extract_coverage('/data3/isaac/early_diff/NEW/quantification/te_count/hg38_rmsk.bed','/data3/isaac/early_diff/NEW/plots/TE_genome_coverage', 3099750718)
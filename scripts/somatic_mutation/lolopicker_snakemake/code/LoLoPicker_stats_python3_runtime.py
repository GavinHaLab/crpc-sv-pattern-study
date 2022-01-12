#!/app/easybuild/software/Python/3.7.4-foss-2016b/bin/python -u

import sys, getopt
#import re
#import numpy
from scipy import stats
#from scipy.stats import binom
#from scipy.cluster.vq import kmeans2, whiten
#from pylab import plot,show
#import pylab
from datetime import datetime

def main(argv):
	basecov = 30000000
	SNPcut = 1

	if len(sys.argv) < 2:
		print('usage: LoLoPicker_stats.py -o <outputpath>')
		sys.exit(1)
	try:
		opts, args = getopt.getopt(argv,"ho:s:g:c:m", ["help","outputpath=","intervalsize=", "genome=", "SNPcutoff=", "method="])
	except getopt.GetoptError:
		print('usage: LoLoPicker_stats.py -o <outputpath>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('usage: LoLoPicker_stats.py -o <outputpath>')
			sys.exit()
		elif opt in ("-o", "--outputpath"):
			outputpath = arg
            ######changed on 12/9/19 by Anna to work with snakemake
			#inputvariant = outputpath + "/control_stats.txt"
			inputvariant = outputpath + "/merged_control_stats.txt"
			statsfile = outputpath+"/stats_calls.txt"
			rejectfile = outputpath+"/reject_calls.txt"
		elif opt in ("-s", "--intervalsize"):
			basecov = arg
	#	elif opt in ("-g", "--genome"):
	#		basecov = 3000000000
		elif opt in ("-c", "--SNPcutoff"):
			SNPcut = arg
		elif opt in ("-m", "--method"):
			correction = arg

	return inputvariant, statsfile, rejectfile, correction, int(basecov), int(SNPcut)

if __name__ == '__main__':
	(inputvariant, statsfile, rejectfile, correction, basecov, SNPcut) = main(sys.argv[1:])

	stats_hash = {}
	varfile = open(inputvariant)
    #Commented out by Anna on 11/25/19 so that we can have headerless inputs,
    #which helps our snakemake
	#next(varfile)
	p_vals=[]

	ftemp1 = open(statsfile, 'w')
	ftemp2 = open(rejectfile, 'w')
	for variant in ( raw.strip().split() for raw in varfile ):
		#print(variant[0], variant[1], "at", datetime.now(), flush = True)
		c_alt_total = 0; c_ref_total = 0; SNP = 0
		c_alf = []; c_alf_info = []
		p_value = 1

		(chr, pos, ref, alt, t_refcount, t_alt_count, n_refcount, n_alt_count, judge, control_info) = variant
		k = (chr+"\t"+str(pos)+"\t"+ref+"\t"+alt)
		stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), str(judge), float(p_value)
		if control_info != "NA":
			controls = control_info.split(',')
			for j in controls:
				(c_alt, c_ref, c_sampleID) = j.split(':')
				c_total = int(c_alt) + int(c_ref)
				if c_total < 10:
					c_alt_total = c_alt_total + int(c_alt)
					c_ref_total = c_ref_total + int(c_ref)
				else:
					alf = int(c_alt)/c_total
					if alf == 0:
						c_alt_total = c_alt_total + int(c_alt)
						c_ref_total = c_ref_total + int(c_ref)
					elif alf >= 0.5:
						SNP += 1
						if SNP > SNPcut:
							stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), "possible_SNP", float(p_value)
							break
					else:
						c_alf.append(alf)
						c_alf_info.append([c_alt, c_ref, c_sampleID])

			#perform kmeans clustering
			#if (len(c_alf) > 3) and ('possible_SNP' not in str(stats_hash[k])):
			#	while True:
			#		km_results = kmeans2(c_alf, 2, iter=10, thresh=0)
			#		km = km_results[0]; alfinfo =str(km_results[1])
			#		cluster1 = km[0]; cluster2 = km[1]
			#		if cluster1 > cluster2:
			#			if (cluster1 < 1) and (cluster2 > 0):
			#				break
			#		if cluster1 < cluster2:
			#			if (cluster2 < 1) and (cluster1 > 0):
			#				break
			#	if (cluster1 > 0.1) and (cluster2 > 0.1):
			#		stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), "possible_SNP", float(p_value)
			#	else:
			#			#no obvious k clusters in the population, only use alf < 0.1 for binom.
			#		if abs(cluster1 - cluster2) < 0.2:
			#			l = 0
			#			for alf in c_alf:
			#				if alf < 0.1:
			#					(c_alt, c_ref, c_sampleID) = c_alf_info[l]
			#					c_alt_total = c_alt_total + int(c_alt)
						#			  				c_ref_total = c_ref_total + int(c_ref)
			#				l += 1
			#		else:
			#				#with two clusters, use lager cluster to filter SNP and use smaller group for binom
			#			if cluster1 > cluster2:
			#				if alfinfo.count('1') > SNPcut:
			#					stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), "possible_SNP", float(p_value)
			#				else:
			#					l = 0
			#					for clust in alfinfo[1::2]:
			#						if clust == '1':
						#									 	 		(c_alt, c_ref, c_sampleID) = c_alf_info[l]
			#							c_alt_total = c_alt_total + int(c_alt)
			#							c_ref_total = c_ref_total + int(c_ref)
			#							l += 1
			#			elif cluster2 > cluster1:
			#				if alfinfo.count('0') > 3:
			#					stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), "possible_SNP", float(p_value)
			#				else:
			#					l = 0
			#					for clust in alfinfo[1::2]:
			#						if clust == '0':
			#							(c_alt, c_ref, c_sampleID) = c_alf_info[l]
			#							c_alt_total = c_alt_total + int(c_alt)
			#							c_ref_total = c_ref_total + int(c_ref)
			#							l += 1

			c_total = int(c_alt_total) + int(c_ref_total)
			t_total = int(t_alt_count) + int(t_refcount)
			n_total = int(n_alt_count) + int(n_refcount)

			if 'possible_SNP' not in str(stats_hash[k]):
				#perform binomial test.##

				if c_total == 0:
					c_alf_all = 1/1000
				else:
				 	c_alf_all = int(c_alt_total)/int(c_total)
				#	if c_alf_all < 1/1000:
				#		c_alf_all = 1/1000

				if c_alf_all == 0:
					c_alf_all = 1/1000
				if int(t_alt_count) >= 150:
					p_value = 0
				else:
					pro = 0
					critical_val = 0
					p_value = stats.binom_test(t_alt_count, t_total, c_alf_all, 'greater')
					for i in range(0, int(t_total)+int(n_total)):
						pro = pro + stats.binom.pmf(i, int(t_total)+int(n_total), c_alf_all)
						if pro > 0.95:
							critical_val = i
							break
					t_power = float(stats.binom_test(critical_val, int(t_total)+int(n_total), 0.3, 'greater'))
					n_power = float(stats.binom_test(critical_val, int(t_total)+int(n_total), 0.5, 'greater'))

					if t_power < 0.8 or n_power < 0.8:
						stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), str(c_alt_total)+":"+str(c_total)+":uncovered", 1
					else:
						stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), str(c_alt_total)+":"+str(c_total), float(p_value)
					#stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), str(c_alt_total)+":"+str(c_total)+":NA", float(p_value)
			else:
				stats_hash[k] = str(t_refcount), str(t_alt_count), str(n_refcount), str(n_alt_count), str(c_alt_total)+":"+str(c_total)+":possibleSNP", 1

	stats_sorted = sorted(list(stats_hash.items()), key=lambda kv: kv[1][5], reverse=True)
	j=1; pre_p = 1

	print("#chr\tpos\tref\talt\ttumor_ref_reads\ttumor_alt_reads\ttumor_alf\treference_ref_reads\treference_alt_reads\tcontrol_info(Readcount_alt:Readcount_ref)\tp_value\tcorrected_p\tjudgement", file=ftemp1)
	print("#chr\tpos\tref\talt\ttumor_ref_reads\ttumor_alt_reads\ttumor_alf\treference_ref_reads\treference_alt_reads\tinfo\tp_value\tcorrected_p\tjudgement", file=ftemp2)

	if correction == "Bonferroni":

		for var in stats_sorted:
			if '#' not in var[0]:
				p = var[1][-1]
				q = float(p) * int(basecov)
				tot = int(var[1][0])+int(var[1][1])
				if q > 1:
					q = 1

				if tot!=0:
					flt = int(var[1][1])/tot
					if q < 0.05:
						results=(str(var[0])+"\t"+str(var[1][0])+"\t"+str(var[1][1])+"\t"+str(flt)+"\t"+str(var[1][2])+"\t"+str(var[1][3])+"\t"+str(var[1][4])+"\t"+str(p)+"\t"+str(q)+"\t"+"real")
						print(results, file=ftemp1)
					else:
						results=(str(var[0])+"\t"+str(var[1][0])+"\t"+str(var[1][1])+"\t"+str(flt)+"\t"+str(var[1][2])+"\t"+str(var[1][3])+"\t"+str(var[1][4])+"\t"+str(p)+"\t"+str(q)+"\t"+"reject")
						print(results, file=ftemp2)



	if correction == "FDR":

		lastq = 1
		i = 1

		for var in stats_sorted:
			if '#' not in var[0]:
				p = var[1][-1]
				q = float(p) * int(basecov)/i
				tot = int(var[1][0])+int(var[1][1])
				if q > lastq:
					q = lastq

				lastq = q
				i = i+1

				if tot!=0:
					flt = int(var[1][1])/tot
					if q < 0.001:
						results=(str(var[0])+"\t"+str(var[1][0])+"\t"+str(var[1][1])+"\t"+str(flt)+"\t"+str(var[1][2])+"\t"+str(var[1][3])+"\t"+str(var[1][4])+"\t"+str(p)+"\t"+str(q)+"\t"+"real")
						print(results, file=ftemp1)

					else:
						results=(str(var[0])+"\t"+str(var[1][0])+"\t"+str(var[1][1])+"\t"+str(flt)+"\t"+str(var[1][2])+"\t"+str(var[1][3])+"\t"+str(var[1][4])+"\t"+str(p)+"\t"+str(q)+"\t"+"reject")
						print(results, file=ftemp2)




	ftemp1.close()
	ftemp2.close()
	print("Done LoLoPicker!", flush = True)

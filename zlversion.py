import pybedtools
from pyfasta import Fasta
import sys
from bx.intervals.intersection import Intersecter, Interval

##########################################################################################
################## MAKE SURE THE ORDER FOR ALL OF SPECIES SHOULD BE SAME! ################
##########################################################################################
zmgff = pybedtools.BedTool("cold2_gff3/Zm284v3.gff3")
#osgff = pybedtools.BedTool("xlai_gff3/Osativa.gff3")
sigff = pybedtools.BedTool("cold2_gff3/Si312.gff3")
sbgff = pybedtools.BedTool("cold2_gff3/Sb313.gff3")
#zmgff = pybedtools.BedTool("GFF_files/")
zmfasta = Fasta("cold2_fasta/Zm284v3.fa")
#osfasta = Fasta("xlai_fasta/Osativa.fa")
sifasta = Fasta("cold2_fasta/Si312.fa")
sbfasta = Fasta("cold2_fasta/Sb313.fa")
#osfasta = Fasta("Fasta_files/Oryza_sativa_v6.1_dna.fasta")

fastas = [zmfasta,sifasta,sbfasta]
gffs = [zmgff,sigff,sbgff]

a = sys.argv[1]
fh = open(a,'r')
fh.readline()
sdict = []
#syn_os = set([])
syn_zm = set([])
syn_si = set([])
syn_sb = set([])
for line in fh:
	new = line.strip().split(',')
	if new[2].startswith('N'):continue
#	if new[1].startswith('N'):continue
	if new[3].startswith('N'):continue
	if new[0].startswith('N'):continue
#	syn_zm.add(new[1])	
	syn_zm.add(new[2])
	syn_si.add(new[3])
	syn_sb.add(new[0])
syn = [syn_zm,syn_si,syn_sb]

syn_region = {}
mdict = {}
for i in gffs:
	k = gffs.index(i)
	syn_region[k] = {}
	for y in i:
		if k == 0:
			if y.name.startswith('GRMZM'):
				y.name = y.name.split('.')[0]
			else:
				y.name = '_'.join(y.name.split('.')[:2])
		if k == 1 or k == 2:
			y.name = '.'.join(y.name.split('.')[:2])
		if y[2] != 'gene':continue
		if y.name not in syn[k]:continue
		if y[0] not in syn_region[k]:
			syn_region[k][y[0]] = Intersecter()	
		syn_region[k][y[0]].add_interval(Interval(int(y[3]),int(y[4]),y.name))
		if y.name not in mdict:
			mdict[y.name] = []
#		if '_' in y.name:continue
#		if 'scaffold' in y.chrom:continue
		mdict[y.name].append(y.chrom)
		mdict[y.name].append(int(y[3]))
		mdict[y.name].append(int(y[4]))
		mdict[y.name].append(y.strand)
		
cand = open(sys.argv[2])
count = 0
for line in cand:
	count += 1
	new = line.strip().split(' ')
	myanns = []
	myseqs = []
	string = []
	for i in new:
		mygene = i 
		mychr = mdict[i][0]
		mystart = mdict[i][1]
		mystop = mdict[i][2]
		mystrand = mdict[i][3]
		k = new.index(i)
		a = syn_region[k][mychr].find(mystart-10000,mystart-1)
		b = syn_region[k][mychr].find(mystop+1,mystop+10000)
		st = []
		sp = []
		if len(a) > 0:
			for c in a:
			#	print c.value
				sp.append(c.end)
			m = max(sp)+1
			myst = m
		else:
			myst = mystart-10000
			if myst < 0:
				myst = 1
		if len(b) > 0:
			for c in b:
			#	print c.value
				st.append(c.start)
			n = min(st)-1
			mysp = n
		else:
			mysp = mystop+10000
		string = [mygene,mychr,mystrand,str(myst),str(mystart),str(mystop),str(mysp)]
		myanns.append("> {0}".format(' '.join(string)))
		myseqs.append(">{0}\n{1}".format(mygene,fastas[k].sequence({"chr":mychr,'start':myst,'stop':mysp,'strand':mystrand})))
	out = open("outputs/gene_extract_{0}.fasta".format(str(count).zfill(5)),'w')
	for q in myanns:
		out.write(q+'\n')
	for d in myseqs:
		out.write(d+'\n')
	out.close()
fh.close()
cand.close()

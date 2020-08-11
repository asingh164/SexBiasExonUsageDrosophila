# Python script originally written by Tyler Kent 2017
# Modified by Amardeep Singh (Modifications noted)

#from collections import OrderedDict
import argparse
import pandas as pd

def arguments():
	parser  = argparse.ArgumentParser(description="Gets 4fold and 0fold sites from bedtools getfasta of CDS annotations")
	parser.add_argument("-i", "--input", help="path for input", required=True)
	parser.add_argument("-o","--output", help="path prefix for output, e.g. ../Data/thaliana", required=True)
	args = parser.parse_args()
	return(args)

def order_check(lines,outpath,Degenerate):

	'''
	Loop over lines and check if next line's start position
	overlaps with previous line's span.
	Send overlapping sections or single nonoverlapping lines
	to function to choose the longest line to write to file.
	'''

#	overlaps = OrderedDict()
#	overlaps = pd.DataFrame(columns = ['line','start','end','span'])
	last = lines[-1]
	fniter = 0
	iter = 0

	for line in lines:
		#set intial line to previous line info for the first comparison
		if iter==0:
			sline = line.split()[0].split(":")
			pchr = sline[2]
			pstart = sline[3].split("-")[0]
			pend = sline[3].split("-")[1].split("(")[0]
			pspan = int(pend)-int(pstart)
			psense = sline[3].split("(")[1].rstrip(")")
			pname = sline[0]
			pseq = line.split()[1]
			overlaps = pd.DataFrame(columns = ['seq','chr','start','end','span','sense','name'])
			overlaps = overlaps.append(pd.Series([pseq,pchr,pstart,pend,pspan,psense,pname],index=['seq','chr','start','end','span','sense','name']),ignore_index=True)
#			check=False
#			newkey = "-".join(pstart,pend,pspan)
#			overlaps[newkey]  = line

		else:
			chr = line.split()[0].split(":")[2]

			#only split and compare things if they're on the same chromosome
			if chr == pchr:
				sline = line.split()[0].split(":")
				start = sline[3].split("-")[0]
				end = sline[3].split("-")[1].split("(")[0]

				if int(start) < int(pend):
					span = int(end)-int(start)
					sense = sline[3].split("(")[1].rstrip(")")
					name = sline[0]
					seq = line.split()[1]
					lineseries = pd.Series([seq,chr,start,end,span,sense,name],index=['seq','chr','start','end','span','sense','name'])
					overlaps = overlaps.append(lineseries,ignore_index=True)
					pstart = start
					pend = end
					pspan = span
					pchr = chr
#					overlaps[newkey] = line

					if line == last:
						codons(Degenerate,overlaps,fniter,outpath)
						#keep_longest(overlaps,outpath,fniter)
						fniter+=1
#						check=True

				else:
					codons(Degenerate,overlaps,fniter,outpath)
					#keep_longest(overlaps,outpath,fniter)
					fniter+=1
#					check=True

					sline = line.split()[0].split(":")
					pchr = sline[2]
					pstart = sline[3].split("-")[0]
					pend = sline[3].split("-")[1].split("(")[0]
					pspan = int(pend)-int(start)
					psense = sline[3].split("(")[1].rstrip(")")
					pname = sline[0]
					pseq = line.split()[1]
					overlaps = pd.DataFrame(columns = ['seq','chr','start','end','span','sense','name'])
					overlaps = overlaps.append(pd.Series([pseq,pchr,pstart,pend,pspan,psense,pname],index=['seq','chr','start','end','span','sense','name']),ignore_index=True)
#					lineseries = pd.Series([pseq,pchr,pstart,pend,pspan,psense,pname],index=['seq','chr','start','end','span','sense','name'])
#					overlaps = overlaps.append(lineseries,ignore_index=True)

					if line == last:
						codons(Degenerate,overlaps,fniter,outpath)
						fniter+=1
#						outfile = open(outpath,	'a')
#						outfile.write(line)
#						outfile.close()

			else:
				codons(Degenerate,overlaps,fniter,outpath)
				#keep_longest(overlaps,outpath,fniter)
				fniter+=1
#				check=True

				sline = line.split()[0].split(":")
				pchr = sline[2]
				pstart = sline[3].split("-")[0]
				pend = sline[3].split("-")[1].split("(")[0]
				pspan = int(pend)-int(start)
				psense = sline[3].split("(")[1].rstrip(")")
				pname =	sline[0]
				pseq = line.split()[1]
				overlaps = pd.DataFrame(columns = ['seq','chr','start','end','span','sense','name'])
				overlaps = overlaps.append(pd.Series([pseq,pchr,pstart,pend,pspan,psense,pname],index=['seq','chr','start','end','span','sense','name']),ignore_index=True)
#				lineseries = pd.Series([pseq,pchr,pstart,pend,pspan,psense,pname],index=['seq','chr','start','end','span','sense','name'])
#				overlaps = overlaps.append(lineseries,ignore_index=True)


				if line == last:
						codons(Degenerate,overlaps,fniter,outpath)
						fniter+=1
#						outfile = open(outpath,	'a')
#						outfile.write(line)
#						outfile.close()

		iter+=1


def keep_longest(overlaps,outpath,iter):

	'''
	Take dictionary of span lengths, get the longest one,
	print the corresponding line to the output file.
	Append to file unless the iteration is the first.
	'''
	overlaps['span'] = overlaps['span'].astype(float)
	df = overlaps.sort_values(['span'], ascending=False)
	line = df['line'].iloc[0]
#	keys = list(overlaps.keys())
#	intkeys = list(map(int, keys))
#	keylist = [x.split("-")[2]

#	longest = max(intkeys)

#	longestline = overlaps[str(longest)]

	if iter == 0:
		outfile = open(outpath, 'w')
		outfile.write(line)

	else:
		outfile = open(outpath, 'a')
		outfile.write(line)

	outfile.close()

def endpos(x):

	'''
	Redefine end position if pos sense not divisible by 3
	'''

	if (int(x['end'])-int(x['start'])) % 3 == 0 and x['sense'] == "+":
		return(int(x['end']))

	elif (int(x['end'])-int(x['start'])) % 3 == 1 and x['sense'] == "+":
		return(int(x['end'])+2)

	elif (int(x['end'])-int(x['start'])) % 3 == 2 and x['sense'] == "+":
		return(int(x['end'])+1)

	elif x['sense'] == "-":
		return(int(x['end']))


def startpos(x):

	'''
	Redefine start position if neg sense not divisible by 3
	'''

	if (int(x['end'])-int(x['start'])) % 3 == 0 and x['sense'] == "-":
		return(int(x['start']))

	elif (int(x['end'])-int(x['start'])) % 3 == 1 and x['sense'] == "-":
		return(int(x['start'])-2)

	elif (int(x['end'])-int(x['start'])) % 3 == 2 and x['sense'] == "-":
		return(int(x['start'])-1)

	elif x['sense'] == "+":
		return(int(x['start']))

def positions(x):

	'''
	Get lists of positions from seq data dependent on sense
	'''

	newstart = x['start'] + 2
	newend = x['end'] + 1

	if x['sense'] == "-":
		test = list(range(newstart,newend,3))[::-1]
		test = [x - 2 for x in test]

	elif x['sense'] == "+":
		test = list(range(newstart,newend,3))

	return(test)

def codons(Degenerate, lines, fniter, outpath):

	'''
	Split CDS into codons. Get lists of
	4fold and 0fold sites.
	Send each to another fn to drop
	non-overlapping site calls in
	overlapping CDS annotations.
	'''

	#split seq into codons
	seq = pd.Series(lines['seq'].apply(lambda x: list(x[i:i+3].upper() for i in range(0, len(x), 3))))

	#get 3rd codon positions
	possense = lines[['sense','start','end']]
	possense['end'] = possense.apply(endpos,axis=1)
	possense['start'] = possense.apply(startpos,axis=1)
#	thirdposseries = possense.apply(positions,axis=1)
#	print(thirdposseries.__class__)
	thirdposseries=pd.Series()
	for i in range(len(possense.index)):
		l = positions(possense.iloc[i])
		thirdposseries=thirdposseries.append(pd.Series([l]))
	thirdposseries=thirdposseries.reset_index(drop=True)

	#get 1st codon positions
	firstposseries = thirdposseries.apply(lambda x: [n-2 for n in x])

	#get 2nd codon positions
	secondposseries = thirdposseries.apply(lambda x: [n-1 for n in x])

	#combine data into single df
	senselist = lines['sense'].tolist()
	startlist = lines['start'].tolist()
	endlist = lines['end'].tolist()
	namelist = lines['name'].tolist()
	chr = lines['chr'][0]

	df = pd.DataFrame(columns = ['chr','codons','first','second','third','sense','start','end','name'])
	df = df.append(pd.Series([chr,seq,firstposseries,secondposseries,thirdposseries,senselist,startlist,endlist,namelist], index = ['chr','codons','first','second','third','sense','start','end','name']), ignore_index = True)

	pull_degenerates(Degenerate,df,fniter,outpath)

def pull_degenerates(Degenerate,df,fniter,outpath):

	'''
	Get full lists of degenerate codons
	from each overlapping line.
	'''

	#loop over overlapping lines

	iter = 0
	for i in range(len(df['codons'][0])):
		seqstacked = pd.Series(df['codons'][0][i])
		firststacked = pd.Series(df['first'][0][i])
		secondstacked = pd.Series(df['second'][0][i])
		thirdstacked = pd.Series(df['third'][0][i])

		dat = pd.concat([seqstacked,firststacked,secondstacked,thirdstacked],axis=1)
		dat.columns = ['codon','first','second','third']

		if iter == 0:
			fold4_p = dat[dat['codon'].isin(Degenerate.get('4fold'))]
			fold4_p = fold4_p[['third']].dropna().reset_index(drop=True)
			fold4_p.columns = ['pos']

			fold0_1p = dat[dat['codon'].isin(Degenerate.get('0fold_1'))]
			fold0_1p = fold0_1p[['first']].dropna().reset_index(drop=True)
			fold0_1p.columns=['pos']

			fold0_2p = dat[dat['codon'].isin(Degenerate.get('0fold_2'))]
			fold0_2p = fold0_2p[['second']].dropna().reset_index(drop=True)
			fold0_2p.columns=['pos']

			fold0_12p = dat[dat['codon'].isin(Degenerate.get('0fold_12'))]
			fold0_12_first = pd.Series(fold0_12p['first'])
			fold0_12_second = pd.Series(fold0_12p['second'])
			fold0_12p = pd.DataFrame(pd.concat([fold0_12_first,fold0_12_second],axis=0).dropna().reset_index(drop=True))
			fold0_12p.columns=['pos']

			fold0_allp = dat[dat['codon'].isin(Degenerate.get('0fold_all'))]
			fold0_allp = pd.DataFrame(pd.concat([pd.Series(fold0_allp['first']),pd.Series(fold0_allp['second']),pd.Series(fold0_allp['third'])],axis=0).dropna().reset_index(drop=True))
			fold0_allp.columns=['pos']

			fold0_p = pd.concat([fold0_1p,fold0_2p,fold0_12p,fold0_allp],axis=0).reset_index(drop=True)

			pstart = int(df['start'][0][0])
			pend = int(df['end'][0][0])
			chr = str(df['chr'][0])
			name = str(df['name'][0][0]).split('.')[0]

		else:
			fold4 = dat[dat['codon'].isin(Degenerate.get('4fold'))]
			fold4 = fold4[['third']].dropna().reset_index(drop=True)
			fold4.columns = ['pos']

			fold0_1 = dat[dat['codon'].isin(Degenerate.get('0fold_1'))]
			fold0_1 = fold0_1[['first']].dropna().reset_index(drop=True)
			fold0_1.columns = ['pos']

			fold0_2 = dat[dat['codon'].isin(Degenerate.get('0fold_2'))]
			fold0_2 = fold0_2[['second']].dropna().reset_index(drop=True)
			fold0_2.columns=['pos']

			fold0_12 = dat[dat['codon'].isin(Degenerate.get('0fold_12'))]
			fold0_12_first = pd.Series(fold0_12['first'])
			fold0_12_second = pd.Series(fold0_12['second'])
			fold0_12 = pd.DataFrame(pd.concat([fold0_12_first,fold0_12_second],axis=0).dropna().reset_index(drop=True))
			fold0_12.columns=['pos']

			fold0_all = dat[dat['codon'].isin(Degenerate.get('0fold_all'))]
			fold0_all = pd.DataFrame(pd.concat([pd.Series(fold0_all['first']),pd.Series(fold0_all['second']),pd.Series(fold0_all['third'])],axis=0).dropna().reset_index(drop=True))
			fold0_all.columns=['pos']

			fold0 = pd.concat([fold0_1,fold0_2,fold0_12,fold0_all],axis=0).reset_index(drop=True)
			fold0.columns = ['pos']

			if int(df['start'][0][i]) < pstart:
				lefthang0 = fold0[(fold0['pos'] >= int(df['start'][0][i])) & (fold0['pos'] < pstart)]
				lefthang4 = fold4[(fold4['pos'] >= int(df['start'][0][i])) & (fold4['pos'] < pstart)]

			elif int(df['start'][0][i]) > pstart:
				lefthang0 = fold0_p[(fold0_p['pos'] >= pstart) & (fold0_p['pos'] < int(df['start'][0][i]))]
				lefthang4 = fold4_p[(fold4_p['pos'] >= pstart) & (fold4_p['pos'] < int(df['start'][0][i]))]

			else:
				lefthang0 = pd.DataFrame(columns=['pos'])
				lefthang4 = pd.DataFrame(columns=['pos'])

			if int(df['end'][0][i]) < pend:
				righthang0 = fold0[(fold0['pos'] > int(df['start'][0][i])) & (fold0['pos'] <= pend)]
				righthang4 = fold4[(fold4['pos'] > int(df['start'][0][i])) & (fold4['pos'] <= pend)]

			elif int(df['end'][0][i]) > pend:
				righthang0 = fold0[(fold0['pos'] > pstart) & (fold0['pos'] <= int(df['start'][0][i]))]
				righthang4 = fold4[(fold4['pos'] > pstart) & (fold4['pos'] <= int(df['start'][0][i]))]

			else:
				righthang0 = pd.DataFrame(columns=['pos'])
				righthang4 = pd.DataFrame(columns=['pos'])

			#keep overlapping sites between overlapping CDS
			inner4 = fold4_p[fold4_p['pos'].isin(fold4['pos'])].dropna()
			inner0 = fold0_p[fold0_p['pos'].isin(fold0['pos'])].dropna()

			#concat with left and right overhangs and reset the 'p' dataframes
			fold4_p = pd.concat([lefthang4,inner4,righthang4])
			fold0_p = pd.concat([lefthang0,inner0,righthang0])

			if int(df['start'][0][i]) < pstart:
				pstart=int(df['start'][0][i])
			if int(df['end'][0][i]) > pend:
				pend=int(df['end'][0][i])

		iter+=1

	#drop 4folds that line up with 0folds on opposite strand
	fold4_p = fold4_p[~fold4_p['pos'].isin(fold0_p['pos'])].dropna()

	#sort both
	fold4_p = fold4_p.sort_values(by=['pos']).reset_index(drop=True)
	fold0_p = fold0_p.sort_values(by=['pos']).reset_index(drop=True)

	#make vector of chr, name
	chr4list = pd.Series([chr for i in range(len(fold4_p))])
	name4list = pd.Series([name for i in range(len(fold4_p))])
	chr0list = pd.Series([chr for i in range(len(fold0_p))])
	name0list = pd.Series([name for i in range(len(fold0_p))])

	#combine
	fold4_p_end = pd.Series(fold4_p['pos'] +1)
	fold0_p_end = pd.Series(fold0_p['pos'] +1)

	fold4_bed = pd.concat([chr4list,fold4_p['pos'],fold4_p_end,name4list],axis=1)
	fold0_bed = pd.concat([chr0list,fold0_p['pos'],fold0_p_end,name0list],axis=1)

	fold4_bed.columns=['chr','start','end','name']
	fold0_bed.columns=['chr','start','end','name']
	fold4_bed[['start','end']] = fold4_bed[['start','end']].astype(int)
	fold0_bed[['start','end']] = fold0_bed[['start','end']].astype(int)

	if fniter==0:
		fold4_bed.to_csv(outpath+'_4fold.bed',header=False,index=False,sep='\t')
		fold0_bed.to_csv(outpath+'_0fold.bed',header=False,index=False,sep='\t')

	else:
		fold4_bed.to_csv(outpath+'_4fold.bed',header=False,index=False,mode='a',sep='\t')
		fold0_bed.to_csv(outpath+'_0fold.bed',header=False,index=False,mode='a',sep='\t')


def define_sites():

	'''
	Create look-up dictionary of degenerate codons
	'''

	Degenerate = {
		"4fold":{"TCT","TCC","TCG","TCA","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
			"CGT","CGC","CGA","CGG","ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG",
			"GCT","GCC","GCA","GCG","GGT","GGC","GGA","GGG"
		},

		"0fold_1":{"TGA","TAG","TAA","TCA"
		},

		"0fold_2":{"AGG","AGA","CGA","CGG","CGC","CGT","CTG","CTC","CTA","CTT","TTG","TTA"
		},

		"0fold_12":{"GGA","GGC","GGG","GGT","GAG","GAA","GAT","GAC","GCA","GCC","GCA",
			"GCT","GTA","GTG","GTC","GTT","AAG","AAA","AAC","AAT","ACG","ACC","ACT",
			"ACA","ATA","ATC","ATT","CAG","CAA","CAT","CAC","CCG","CCC","CCA","CCT",
			"TGC","TGT","TAC","TAT","TTC","TTT","TCT","TCC","TCG","AGC","AGT","GCG"

		},

		"0fold_all":{"ATG","TGG"
		}
	}

	return(Degenerate)

#####################
# BEGIN
#####################

args = arguments()

inpath = args.input
outpath = args.output

infile = open(inpath, 'r')
lines = infile.readlines()

#get degenerate codons
Degenerate = define_sites()

#pass all lines to function
order_check(lines,outpath,Degenerate)

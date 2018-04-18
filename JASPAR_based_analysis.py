### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
	python JASPAR_based_analysis.py
	--jaspar <FULL_PATH_TO_JASPAR_BED_FILE_WITH_MOTIFS_OF_INTEREST>
	--promoters <FULL_PATH_TO_GFF3_OUTPUT_FILE_FOR_PROMOTER_POSITIONS>
	--results <FULL_PATH_TO_OUTPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import re, sys
from operator import itemgetter
import matplotlib.pyplot as plt

# --- end of imports --- #

def load_CRE_positions( jaspar_file ):
	"""! @brief load all CRE positions from JASPAR file """
	
	raw_CREs = []
	chromosomes = []
	
	with open( jaspar_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			raw_CREs.append( { 'seq_id': parts[0], 'start': int( parts[1] ), 'end': int( parts[2] ), 'id': parts[3], 'score': int( parts[4] ), 'orientation': parts[5] } )
			chromosomes.append( parts[0] )
			line = f.readline()
	sorted_chromosomes = sorted( list( set( chromosomes ) ) )	
	
	CREs = {}
	counter = 0
	for chromosome in sorted_chromosomes:
		values = []
		for CRE in raw_CREs:
			if CRE['seq_id'] == chromosome:
				values.append( CRE )
		CREs.update( { chromosome: values } )
		counter += len( values )
	
	print jaspar_file
	print "number of elements: " + str( counter )
	
	return CREs


def load_promoter_regions( promoter_pos_file ):
	"""! @brief load promoter positions again """
	
	chromosome_name_mapping_table = { 	'NC_000001.11': "chr1",
																		'NC_000002.12': "chr2",
																		'NC_000003.12': "chr3",
																		'NC_000004.12': "chr4",
																		'NC_000005.10': "chr5",
																		'NC_000006.12': "chr6",
																		'NC_000007.14': "chr7",
																		'NC_000008.11': "chr8",
																		'NC_000009.12': "chr9",
																		'NC_000010.11': "chr10",
																		'NC_000011.10': "chr11",
																		'NC_000012.12': "chr12",
																		'NC_000013.11': "chr13",
																		'NC_000014.9': "chr14",
																		'NC_000015.10': "chr15",
																		'NC_000016.10': "chr16",
																		'NC_000017.11': "chr17",
																		'NC_000018.10': "chr18",
																		'NC_000019.10': "chr19",
																		'NC_000020.11': "chr20",
																		'NC_000021.9': "chr21",
																		'NC_000022.11': "chr22",
																		'NC_000023.11': "chrX",
																		'NC_000024.10': "chrY"
																	}
	
	promoter_positions = []
	counter = 0
	with open( promoter_pos_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					promoter_positions.append( { 'id': parts[-1], 'seq_id': chromosome_name_mapping_table[ parts[0] ], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } )
				except KeyError:
					counter += 1
			line = f.readline()
	print "number of missed genes (located on genomic patches): " + str( counter )
	return promoter_positions


def run_analysis( promoter_regions, CREs, all_data_output_file ):
	"""! @brief run analysis to identify cutoff values """
	
	promoter_with_A_B_C = [ [], [], [] ]
	all_data = []
	CRE_types = { "RELA": 0, "RELB":1, "REL":2 }
	
	with open( all_data_output_file, "w", 0 ) as out:
		out.write( "PromoterID\tCRE\tChromosome\tStart\tEnd\tScore\tOrientation\n" )
		for promoter in promoter_regions:
			try:
				candidate_CREs = CREs[ promoter['seq_id'] ]
				for candidate in candidate_CREs:
					if candidate['start'] > promoter['start']:
						if candidate['end'] < promoter['end']:
							promoter_with_A_B_C[ CRE_types[ candidate['id'] ] ].append( promoter['id'] )
							entry = { 	'p_id': promoter['id'],
												'cre': candidate['id'],
												'seq_id': candidate['seq_id'],
												'start': candidate['start'],
												'end': candidate['end'],
												'score': candidate['score'],
												'orientation': candidate['orientation']
										  }
							all_data.append( entry )
							out.write( "\t".join( map( str, [ entry['p_id'], entry['cre'], entry['seq_id'], entry['start'], entry['end'], entry['score'], entry['orientation'] ] ) ) + '\n' )
			except KeyError:
				print "KeyError: " + str( promoter['seq_id']  )
	
	output_files = [ 	all_data_output_file + "_A.txt",
								all_data_output_file + "_B.txt",
								all_data_output_file + "_C.txt" ]
	
	for idx, CRE in enumerate( promoter_with_A_B_C ):
		print len( CRE )
		with open( output_files[ idx ], "w" ) as out:
			out.write( "\n".join( list( set( CRE ) ) ) + '\n' )


def load_results( all_data_output_file ):
	"""! @brief load all data from result file to save time by avoiding time consuming computing of data """
	
	results = []
	
	with open( all_data_output_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			results.append( { 	'p_id': parts[0].strip(),
											'cre': parts[1].strip(),
											'seq_id': parts[2].strip(),
											'start': int( parts[3].strip() ),
											'end': int( parts[4].strip() ),
											'score': int( parts[5].strip()),
											'orientation': parts[6].strip()
										} )
			line = f.readline()
	return results


def quality_score_distributions( data, output_dir ):
	
	rela_scoreA = []
	relb_scoreA = []
	relc_scoreA = []
	
	rela_scoreB = []
	relb_scoreB = []
	relc_scoreB = []
	
	for each in data:
		if each['cre'] == "RELA":
			rela_scoreA.append( each['scoreA'] )
			rela_scoreB.append( each['scoreB'] )
		elif each['cre'] == "RELB":
			relb_scoreA.append( each['scoreA'] )
			relb_scoreB.append( each['scoreB'] )
		else:
			relc_scoreA.append( each['scoreA'] )
			relc_scoreB.append( each['scoreB'] )
	
	# --- scoreA histograms --- #
	rela_score_hist = output_dir + "rela_scoreA_hist.png"
	fig, ax = plt.subplots()
	ax.hist( rela_scoreA, bins=200 )
	ax.set_title( "RELA scoreA distribution" )
	ax.set_xlabel( "scoreA" )
	ax.set_ylabel( "number of hits" )
	fig.savefig( rela_score_hist, dpi=300 )
	
	relb_score_hist = output_dir + "relb_scoreA_hist.png"
	fig, ax = plt.subplots()
	ax.hist( relb_scoreA, bins=200 )
	ax.set_title( "RELB scoreA distribution" )
	ax.set_xlabel( "scoreA" )
	ax.set_ylabel( "number of hits" )
	fig.savefig( relb_score_hist, dpi=300 )
	
	relc_score_hist = output_dir + "relc_scoreA_hist.png"
	fig, ax = plt.subplots()
	ax.hist( relc_scoreA, bins=200 )
	ax.set_title( "c-REL scoreA distribution" )
	ax.set_xlabel( "scoreA" )
	ax.set_ylabel( "number of hits" )
	fig.savefig( relc_score_hist, dpi=300 )
	
	# --- scoreB histograms --- #
	rela_score_hist = output_dir + "rela_scoreB_hist.png"
	fig, ax = plt.subplots()
	ax.hist( rela_scoreB, bins=200 )
	ax.set_title( "RELA scoreB distribution" )
	ax.set_xlabel( "scoreB" )
	ax.set_ylabel( "number of hits" )
	fig.savefig( rela_score_hist, dpi=300 )
	
	relb_score_hist = output_dir + "relb_scoreB_hist.png"
	fig, ax = plt.subplots()
	ax.hist( relb_scoreB, bins=200 )
	ax.set_title( "RELB scoreB distribution" )
	ax.set_xlabel( "scoreB" )
	ax.set_ylabel( "number of hits" )
	fig.savefig( relb_score_hist, dpi=300 )
	
	relc_score_hist = output_dir + "relc_scoreB_hist.png"
	fig, ax = plt.subplots()
	ax.hist( relc_scoreB, bins=200 )
	ax.set_title( "c-REL scoreB distribution" )
	ax.set_xlabel( "scoreB" )
	ax.set_ylabel( "number of hits" )
	fig.savefig( relc_score_hist, dpi=300 )
	
	plt.close("all")


def get_genes_with_top_elements( data, output_dir ):
	
	# --- get entries sorted by CRE --- #
	a_entries = []
	b_entries = []
	c_entries = []
	
	for each in data:
		if each['cre'] == "RELA":
			a_entries.append( each)
		elif each['cre'] == "RELB":
			b_entries.append( each )
		else:
			c_entries.append( each )
	
	scoreA_sorted_A = sorted( a_entries, key=itemgetter('scoreA', 'scoreB') )
	scoreB_sorted_A = sorted( a_entries, key=itemgetter('scoreB', 'scoreA') )
	
	scoreA_sorted_B = sorted( b_entries, key=itemgetter('scoreA', 'scoreB') )
	scoreB_sorted_B = sorted( b_entries, key=itemgetter('scoreB', 'scoreA') )

	scoreA_sorted_C = sorted( c_entries, key=itemgetter('scoreA', 'scoreB') )
	scoreB_sorted_C = sorted( c_entries, key=itemgetter('scoreB', 'scoreA') )
	
	# --- get top CREs (5% FDR) --- #
	a_entries_scoreA = scoreA_sorted_A[ -( int( len( scoreA_sorted_A )*0.05 ) ): ]
	b_entries_scoreA = scoreA_sorted_B[ -( int( len( scoreA_sorted_B )*0.05 ) ): ]
	c_entries_scoreA = scoreA_sorted_C[ -( int( len( scoreA_sorted_C )*0.05 ) ): ]
	
	
	a_entries_scoreB = scoreB_sorted_A[ -( int( len( scoreB_sorted_A )*0.05 ) ): ]
	b_entries_scoreB = scoreB_sorted_B[ -( int( len( scoreB_sorted_B )*0.05 ) ): ]
	c_entries_scoreB = scoreB_sorted_C[ -( int( len( scoreB_sorted_C )*0.05 ) ): ]
	
	
	# --- get overlaps between scoreA and scoreB --- #
	a_file = output_dir + "A_0.05FDR.txt"
	b_file = output_dir + "B_0.05FDR.txt"
	c_file = output_dir + "C_0.05FDR.txt"
	
	
	with open( a_file, "w" ) as out:
		for each in a_entries_scoreA:
			for entry in a_entries_scoreB:
				if each['p_id'] == entry['p_id']:
					out.write( each['p_id'] + '\n' )
	
	
	with open( b_file, "w" ) as out:
		for each in b_entries_scoreA:
			for entry in b_entries_scoreB:
				if each['p_id'] == entry['p_id']:
					out.write( each['p_id'] + '\n' )
	
	
	with open( c_file, "w" ) as out:
		for each in c_entries_scoreA:
			for entry in c_entries_scoreB:
				if each['p_id'] == entry['p_id']:
					out.write( each['p_id'] + '\n' )


def promoter_centered_analysis( results, output_dir ):
	#promoter centered analysis: calculate score per promoter based on all predicted CREs in it
	
	# --- scoreA-based analysis --- #
	a_promoters = {}
	b_promoters = {}
	c_promoters = {}
	
	for each in results:
		if each['cre'] == "RELA":
			try:
				value = a_promoters[ each['p_id'] ]
				value += each['score']
				del a_promoters[ each['p_id'] ]
				a_promoters.update( { each['p_id']: value } )
			except KeyError:
				a_promoters.update( { each['p_id']: each['score'] } )
		
		elif each['cre'] == "RELB":
			try:
				value = b_promoters[ each['p_id'] ]
				value += each['score']
				del b_promoters[ each['p_id'] ]
				b_promoters.update( { each['p_id']: value } )
			except KeyError:
				b_promoters.update( { each['p_id']: each['score'] } )
		
		else:
			try:
				value = c_promoters[ each['p_id'] ]
				value += each['score']
				del c_promoters[ each['p_id'] ]
				c_promoters.update( { each['p_id']: value } )
			except KeyError:
				c_promoters.update( { each['p_id']: each['score'] } )

	a_promoters2 = []
	for key in a_promoters.keys():
		a_promoters2.append( { 'id': key, 'value': a_promoters[ key ] } )
	
	b_promoters2 = []
	for key in b_promoters.keys():
		b_promoters2.append( { 'id': key, 'value': b_promoters[ key ] } )
	
	c_promoters2 = []
	for key in c_promoters.keys():
		c_promoters2.append( { 'id': key, 'value': c_promoters[ key ] } )
	
	hist_rela_promoter_scoreA = output_dir + "RELA_hist_promoter_scoreA.png"
	hist_relb_promoter_scoreA = output_dir + "RELB_hist_promoter_scoreA.png"
	hist_relc_promoter_scoreA = output_dir + "cREL_hist_promoter_scoreA.png"
	
	fig, ax = plt.subplots( )
	ax.hist( a_promoters.values(), bins=200 )
	ax.set_title( "RELA promoters score distribution" )
	ax.set_xlabel( "score" )
	ax.set_ylabel( "number of promoters" )
	fig.savefig( hist_rela_promoter_scoreA, dpi=300 )
	
	fig, ax = plt.subplots( )
	ax.hist( a_promoters.values(), bins=200 )
	ax.set_title( "RELB promoters score distribution" )
	ax.set_xlabel( "score" )
	ax.set_ylabel( "number of promoters" )
	fig.savefig( hist_relb_promoter_scoreA, dpi=300 )

	fig, ax = plt.subplots( )
	ax.hist( a_promoters.values(), bins=200 )
	ax.set_title( "cREL promoters score distribution" )
	ax.set_xlabel( "score" )
	ax.set_ylabel( "number of promoters" )
	fig.savefig( hist_relc_promoter_scoreA, dpi=300 )
	
	# --- combine both analysis results --- #
	a_promoters_005_FDR_A = sorted( a_promoters2, key=itemgetter('value') )[ -( int( len( a_promoters2 )*0.05 ) ): ]
	b_promoters_005_FDR_A = sorted( b_promoters2, key=itemgetter('value') )[ -( int( len( b_promoters2 )*0.05 ) ): ]
	c_promoters_005_FDR_A = sorted( c_promoters2, key=itemgetter('value') )[ -( int( len( c_promoters2 )*0.05 ) ): ]
	
	
	a_file = output_dir + "A_promoter_centered_0.05FDR.txt"
	b_file = output_dir + "B_promoter_centered_0.05FDR.txt"
	c_file = output_dir + "C_promoter_centered_0.05FDR.txt"
	
	with open( a_file, "w" ) as out:
		for each in a_promoters_005_FDR_A:
			out.write( each['id'] + '\t' + str( each['value'] ) + '\n' )
			
	with open( b_file, "w" ) as out:
		for each in b_promoters_005_FDR_A:
			out.write( each['id'] + '\t' + str( each['value'] ) + '\n' )
	
	with open( c_file, "w" ) as out:
		for each in c_promoters_005_FDR_A:
			out.write( each['id'] + '\t' + str( each['value'] ) + '\n' )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	jaspar_file = arguments[ arguments.index('--jaspar')+1 ]	#RELs.bed from JASPAR
	promoter_pos_file = arguments[ arguments.index('--promoters')+1 ]	#output file: promoters.gff3
	all_data_output_file = arguments[ arguments.index('--results')+1 ]	#result file: all_data.txt
	
	output_dir = arguments[ arguments.index('--out')+1 ]	#output directory
	
	# --- need to be run once to get all assignments done (can be skipped later) --- #	
	promoter_regions = load_promoter_regions( promoter_pos_file )
	CREs = load_CRE_positions( jaspar_file )
	run_analysis( promoter_regions, CREs, all_data_output_file )
	
	results = load_results( all_data_output_file )
	
	quality_score_distributions( results, output_dir )
	get_genes_with_top_elements( results, output_dir )
	promoter_centered_analysis( results, output_dir )

if __name__ == '__main__':
	
	if '--jaspar' in sys.argv and '--promoters' in sys.argv and '--results' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

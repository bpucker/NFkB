### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python genome_wide_distribution.py
	--gff <FULL_PATH_TO_GFF3_FILE>
	--results <FULL_PATH_TO_RESULT_FILE>
	--out <FULL_PATH_TO_OUTPUT_FILE>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import sys

# --- end of imports --- #


def load_gene_content_per_chromosome( gff3_file ):
	"""! @brief load number of genes per chromosome """
	
	genes_per_chr = {}
	
	with open( gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					value = genes_per_chr[ parts[0] ] + 1
					del genes_per_chr[ parts[0] ]
					genes_per_chr.update( { parts[0]: value } )
				except KeyError:
					genes_per_chr.update( {  parts[0]:  1} )
			line = f.readline()
	return genes_per_chr


def load_number_of_CREs_per_chromosome( result_file ):
	"""! @brief load elements per chromosome """
	
	RELA_per_chr = {}
	RELB_per_chr = {}
	RELC_per_chr = {}
	
	with open( result_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if "RELA" in parts[1]:
				try:
					value = RELA_per_chr[ parts[ 2 ] ] + 1
					del RELA_per_chr[ parts[ 2 ] ] 
					RELA_per_chr.update( { parts[ 2 ]: value } )
				except KeyError:
					RELA_per_chr.update( { parts[ 2 ]: 1 } )
			elif "RELB" in parts[1]:
				try:
					value = RELB_per_chr[ parts[ 2 ] ] + 1
					del RELB_per_chr[ parts[ 2 ] ] 
					RELB_per_chr.update( { parts[ 2 ]: value } )
				except KeyError:
					RELB_per_chr.update( { parts[ 2 ]: 1 } )
			else:
				try:
					value = RELC_per_chr[ parts[ 2 ] ] + 1
					del RELC_per_chr[ parts[ 2 ] ] 
					RELC_per_chr.update( { parts[ 2 ]: value } )
				except KeyError:
					RELC_per_chr.update( { parts[ 2 ]: 1 } )
			line = f.readline()
	return RELA_per_chr, RELB_per_chr, RELC_per_chr


def construct_genome_wide_distribution_figure( genes_per_chr, chromosome_name_mapping_table, RELA_per_chr, RELB_per_chr, RELC_per_chr, output_file ):
	"""! @brief construct genome wide distribution figure for CREs """
	
	# --- invert mapping table for chromosome names --- #
	new_mapping_table = {}
	for value in chromosome_name_mapping_table.values():
		for key in chromosome_name_mapping_table.keys():
			if chromosome_name_mapping_table[ key ] == value:
				new_mapping_table.update( { value: key } )
	
	# --- construct genome wide CRE distribution plot --- #
	plt.close( "all" )
	fig, ax = plt.subplots( figsize=( 20, 5 ) )
	labels = []
	
	chromosomes = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY" ]
	for idx, key in enumerate( chromosomes ):
		A_value = RELA_per_chr[ key ] /  ( genes_per_chr[ new_mapping_table[ key ] ] / 1000.0)
		ax.plot( [ idx, idx ], [ 0, A_value ], color="black" )
		
		B_value = RELB_per_chr[ key ] /  ( genes_per_chr[ new_mapping_table[ key ] ] / 1000.0)
		ax.plot( [ idx+.1, idx+.1 ], [ 0, B_value ], color="red" )
		
		C_value = RELC_per_chr[ key ] /  ( genes_per_chr[ new_mapping_table[ key ] ] / 1000.0)
		ax.plot( [ idx+.2, idx+.2 ], [ 0, C_value ], color="blue" )
		
		labels.append( key )
	
	
	
	ax.set_xlim( 0, len( labels ) )
	start, end = ax.get_xlim()
	ax.xaxis.set_ticks(np.arange(start, end, 1))
	ax.set_xticklabels( labels )
	ax.set_ylabel( "number of CREs per 1,000 genes" )
	
	ax.set_frame_on(False)
	ax.get_xaxis().tick_bottom()
	
	black_patch = mpatches.Patch( color="black", label="RELA" )
	red_patch = mpatches.Patch( color="red", label="RELB" )
	blue_patch = mpatches.Patch( color="blue", label="cREL" )
	ax.legend( handles=[ black_patch, red_patch, blue_patch ], bbox_to_anchor=(0.5,0.9) )
	
	plt.subplots_adjust( left=0.05, right=0.99, top=0.97, bottom=0.05 )
	
	fig.savefig( output_file )


def main( arguments ):
	"""! @brief run everything """
	
	gff3_file = arguments[ arguments.index( '--gff' )+1 ]	#reference gff3 file: GRCh38_latest_genomic.gff
	result_file = arguments[ arguments.index( '--results' )+1 ]	#result file: all_data.txt
	output_file = arguments[ arguments.index( '--out' )+1 ]	#genome-wide figure: genome_wide.png
	
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
	
	genes_per_chr = load_gene_content_per_chromosome( gff3_file )
	RELA_per_chr, RELB_per_chr, RELC_per_chr = load_number_of_CREs_per_chromosome( result_file )
	
	construct_genome_wide_distribution_figure( genes_per_chr, chromosome_name_mapping_table, RELA_per_chr, RELB_per_chr, RELC_per_chr, output_file )


if __name__ == '__main__':
	if '--gff' in sys.argv and '--results' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

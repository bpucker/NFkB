### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
	python construct_CRE_per_promoter_hist.py
	--results <FULL_PATH_TO_RESULT_FILE>
	--prefix <FULL_PATH_TO_OUTPUT_DIRECTORY>
					"""

import matplotlib.pyplot as plt
import sys, os

# --- end of imports --- #

def main( arguments ):
	"""! @brief construct histograms for CRE distribution ofver promoters """
	
	input_file = arguments[ arguments.index( '--results' )+1 ]	#result file: all_data.txt
	prefix = arguments[ arguments.index( '--prefix' )+1 ]	#prefix for figure files
	
	if prefix[-1] != '/':
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	promoters = []
	rela_promoters = []
	relb_promoters = []
	relc_promoters = []

	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			promoters.append( parts[0] )
			if "RELA" in parts[ 1 ]:
				rela_promoters.append( parts[0] )
			elif "RELB" in parts[1]:
				relb_promoters.append( parts[0] )
			elif"REL" in parts[1]:
				relc_promoters.append( parts[0] )
			line = f.readline()


	promoter_values = []
	for each in list( set( promoters ) ):
		promoter_values.append( promoters.count( each ) )

	rela_promoter_values = []
	for each in list( set( rela_promoters ) ):
		rela_promoter_values.append( rela_promoters.count( each ) )
	with open( prefix + "RELA_genes.txt", "w" ) as out:
		out.write( "\n".join( list( set( rela_promoters ) ) ) )

	relb_promoter_values = []
	for each in list( set( relb_promoters ) ):
		relb_promoter_values.append( relb_promoters.count( each ) )
	with open( prefix + "RELB_genes.txt", "w" ) as out:
		out.write( "\n".join( list( set( relb_promoters ) ) ) )
		
	relc_promoter_values = []
	for each in list( set( relc_promoters ) ):
		relc_promoter_values.append( relc_promoters.count( each ) )
	with open( prefix + "RELC_genes.txt", "w" ) as out:
		out.write( "\n".join( list( set( relc_promoters ) ) ) )

	outfile = prefix + "number_of_binding_sites_per_promoter.png"
	plt.close()
	fig, ax = plt.subplots()
	ax.hist( promoter_values, bins=100, color="black", label="all" )
	ax.hist( rela_promoter_values, bins=100, color="red", label="RELA" )
	ax.hist( relb_promoter_values, bins=100, color="green", label="RELB" )
	ax.hist( relc_promoter_values, bins=100, color="blue", label="cREL" )
	ax.legend()
	fig.savefig( outfile, dpi=300 )


if __name__ == '__main__':
	if '--results' in sys.argv and '--prefix' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

# annotation for TAGseq Experiment
#
######################################################################
#
# assigns 3' UTR's of a given length to reads
#
######################################################################

import os
import argparse

class annApp():

	def __init__(self):
		self.verbose = False

	def start(self, Files, Length):

		directory = os.getcwd()+'/' 

		if not os.path.exists(directory+'edited_gtf/'):   #  creates directory for output files 
			os.makedirs(directory+'edited_gtf/')

		all_files = os.listdir(directory+Files+'/') #  directory  containing all unedited gtf, must be in same directory as annotation.py (you know the drill)

		for file_to_edit in all_files:  #  searches through all files in specified directory for files in gtf format
			if ".gff" in file_to_edit or ".gtf" in file_to_edit:
				with open(directory+Files+'/'+file_to_edit, 'r') as file:
					data = file.readlines()
					data.insert(0,'0\t1\t2\t3\t4\t5\t6\t7\t\"dumby\"; Parent= \"dumby\"\n')  #  dumby line (this program was tricky just live with it)
					file.close()

				# sample line of gtf file

				#fo_apii_R3_nrrl38295_contig_1	CodingQuarry_v2.0	CDS	29841	30254	.	-	0	transcript_id "PGN.00001"; Parent= "PGN.00001";
				#		0							1				2		3		4	5	6	7							8		

				bp_to_add = int(Length)
				
				x = 1

				while x < len(data):  
					i = 1
					con = True
					while con == True:
						if x + i < len(data):  #  for each line in gtf file, counts following lines with identical gene ID (counts number of lines to skip and insert UTR)
							split = data[x].split('\t')
							if "three_prime_UTR" in split[2]:
								x += 1
							elif data[x+i].split('\t')[8].split('Parent=')[1].split('\n')[0] == data[x].split('\t')[8].split('Parent=')[1].split('\n')[0]:
								i += 1
							else:
								i -= 1
								split = data[x].split('\t')
								con = False
						else:
							i -= 1
							con = False

					if 'UTR' not in data[x] and 'UTR' not in data[x+i]:  #  prevents additional UTR lines 
						if split[6] == '-':
							
							if int(split[3])-1-bp_to_add <= 0:
								utr_start = str(1)
							else:
								utr_start = str(int(split[3])-1-bp_to_add)

							if data[x-1].split('\t')[8].split('Parent=')[1].split('\n')[0] != data[x].split('\t')[8].split('Parent=')[1].split('\n')[0]:  #  if '-' strand is first of it's gene ID
								data[x-1] += split[0]+'\t'+split[1]+'\t'+"UTR"+'\t'+utr_start+'\t'+str(int(split[3])-1)+'\t'+split[5]+'\t'+split[6]+'\t'+str(0)+'\t'+split[8]
							elif '\tUTR\t' not in data[x-1]:  #  else if '-' strand Parent already has a UTR line
								data[x-1] += split[0]+'\t'+split[1]+'\t'+"UTR"+'\t'+utr_start+'\t'+str(int(split[3])-1)+'\t'+split[5]+'\t'+split[6]+'\t'+str(0)+'\t'+split[8]
							
						elif 'UTR' not in data[x+i]:  #  prevents addition of extra '+' strand UTR lines
							utr_start = int(data[x+i].split('\t')[4])  #  start of UTR for '+' strand
							data[x+i] += split[0]+'\t'+split[1]+'\t'+"UTR"+'\t'+str(utr_start+1)+'\t'+str(utr_start+1+bp_to_add)+'\t'+split[5]+'\t'+split[6]+'\t'+str(0)+'\t'+split[8]
						
					x += i + 1  #  skips by number of following lines with identical gene ID		

				data[0] = data[0].split('\n')[1]+'\n'  #  corrects file format (removes dumby and blank lines)
				while data[0][:1] == '\n':
					data[0] = data[0][2:]


				with open(directory+'edited_gtf/edited_'+str(bp_to_add)+'_'+file_to_edit, 'w') as file:  #  creates output files
					file.writelines(data)
					file.close()

#####################################################################################
# Help Page and Argument Parsing

def parseArgs():
	#
	# Required and Optional Arguments 
	#
	parser = argparse.ArgumentParser('annotation', description='assigns 3\' UTRs of a specified length to reads', add_help=True,
		epilog='For questions or comments, contact Bradley Jenner <bnjenner@ucdavis.edu>')
	parser.add_argument('-g', help='directory where unedited gtf files are located', dest='Files', type=str)
	parser.add_argument('-l', help='desired length of UTR', dest='Length', type=str)
	
	args = parser.parse_args()

	return args

def execute(args):
	app = annApp()
	return app.start(args.Files, args.Length)

def main():
	"""
	main function
	""" 
	args = parseArgs()
	execute(args)
	
if __name__ == '__main__':
	main()

		

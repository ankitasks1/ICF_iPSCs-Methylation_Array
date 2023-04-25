import os
import sys
import argparse
import time
start_time = time.process_time()
filePath = os.getcwd()
print('Current directory: ' + filePath)


parser = argparse.ArgumentParser(prog='Feature information extractor from gencode',description='Extract features from gtf file')
parser.add_argument('--input', help='Name of the input gtf file', type=argparse.FileType('r'))
parser.add_argument('--feature', help='Name of the feature [eg: gene,exon,transcript]')
parser.add_argument('--output', help='BaseName of the output  file', default='_gencode_out.txt')
parser.add_argument('--species', help='Species name eg. human',  default='human')

# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

#get the name of the inout file
print('Given GTF file : ' + args.input.name)

genome_version = args.input.name.split('.')[1]
print('Genome version : '+genome_version)

#if the file already exist as the same name as output file then delete it
for item in os.listdir(filePath):
   if item.endswith(args.feature + '_' + genome_version + '_' + args.species + '_' + args.output
):
       os.remove(os.path.join(filePath, item))

mygtf = args.input.read().strip().split('\n')
#print(mygtf)

print('''
Author: Ankit Verma
Contact: ankitverma9079@gmail.com
Note: This code will work only with gencode gtf file
''')
def get_info(myinputfile):
   for temp1 in myinputfile:
       temp1 = temp1.split('\t')
       #print(temp1[0])
       #Extract the information of the given gtf file
       if temp1[0].startswith('##description:'):
           print('Description of GTF file: ' + temp1[0].replace('##description:', ''))
       if temp1[0].startswith('##provider:'):
           print('Source of GTF file: ' + temp1[0].replace('##provider:', ''))
       if temp1[0].startswith('##format:'):
           print('Original format of GTF file: ' + temp1[0].replace('##format:', ''))



#Extract the feature, i.e. the column third

def get_genesinfo(myinputfile):
   for temp2 in myinputfile[5:]:
       temp2 = temp2.split('\t')
       #print(temp2)
       if temp2[2] == 'gene':
           chr = temp2[0]
           start = temp2[3]
           end = temp2[4]
           strand = temp2[6]
           featuretype = temp2[2]
           ensmble_id = str(temp2[8].split(' ')[1].strip('";'))
           gene_symbol = str(temp2[8].split(' ')[5].strip('";'))
           #print(chr + '\t' + start + '\t' +end + '\t' +strand + '\t' + featuretype + '\t' + ensmble_id+ '\t' + symbol)
           genes_info = ''.join(chr + '\t' + start + '\t' + end + '\t' + strand + '\t' + featuretype + '\t' + ensmble_id + '\t' + gene_symbol)
           #print(genes_info)
           with open(args.feature + '_' + genome_version + '_' + args.species + '_' + args.output, 'a') as mygeneoutfile:
               mygeneoutfile.write(genes_info + '\n')
               mygeneoutfile.close()

def get_transcriptinfo(myinputfile):
   for temp2 in myinputfile[5:]:
       temp2 = temp2.split('\t')
       # print(temp2)
       if temp2[2] == 'transcript':
           chr = temp2[0]
           start = temp2[3]
           end = temp2[4]
           strand = temp2[6]
           featuretype = temp2[2]
           ensmble_id = str(temp2[8].split(' ')[3].strip('";'))
           gene_symbol = str(temp2[8].split(' ')[7].strip('";'))
           transcript_symbol = str(temp2[8].split(' ')[11].strip('";'))
           #print(chr + '\t' + start + '\t' +end + '\t' +strand + '\t' + featuretype + '\t' + ensmble_id+ '\t' + symbol)
           transcript_info = ''.join(chr + '\t' + start + '\t' + end + '\t' + strand + '\t' + featuretype + '\t' + ensmble_id + '\t' + gene_symbol + '\t' + transcript_symbol)
           #print(transcript_info)
           with open(args.feature + '_' + genome_version  + '_' + args.species + '_' + args.output, 'a') as mytoutfile:
               mytoutfile.write(transcript_info + '\n')
               mytoutfile.close()

def get_exoninfo(myinputfile):
   for temp2 in myinputfile[5:]:
       temp2 = temp2.split('\t')
       # print(temp2)
       if temp2[2] == 'exon':
           chr = temp2[0]
           start = temp2[3]
           end = temp2[4]
           strand = temp2[6]
           featuretype = temp2[2]
           ensmble_id = str(temp2[8].split(' ')[3].strip('";'))
           gene_symbol = str(temp2[8].split(' ')[7].strip('";'))
           exon_symbol = str(temp2[8].split(' ')[15].strip('";'))
           exon_numbert = str(temp2[8].split(' ')[12].strip('";'))
           exon_numberN = str(temp2[8].split(' ')[13].strip('";'))
           #print(chr + '\t' + start + '\t' +end + '\t' +strand + '\t' + featuretype + '\t' + ensmble_id+ '\t' + symbol)
           exon_info = ''.join(chr + '\t' + start + '\t' + end + '\t' + strand + '\t' + featuretype + '\t' + ensmble_id + '\t' + gene_symbol + '\t' + exon_numbert + '_' + exon_numberN + '\t' +  exon_symbol)
           #print(exon_info)
           with open(args.feature + '_' + genome_version + '_' + args.species + '_' + args.output, 'a') as myexonoutfile:
               myexonoutfile.write(exon_info + '\n')
               myexonoutfile.close()

       else:
           pass



if args.feature == 'gene':
   get_info(mygtf)
   get_genesinfo(mygtf)
elif args.feature == 'transcript':
   get_info(mygtf)
   get_transcriptinfo(mygtf)
elif args.feature == 'exon':
   get_info(mygtf)
   get_exoninfo(mygtf)
else:
   print('--feature is missing, please provide appropriate option [gene, transcript, exon]')

print('Output file : ' + args.feature + '_' + genome_version + '_' + args.output)
print('    ')
print('Time taken: ' + str(time.process_time() - start_time))


import argparse
import mass_spectrom_function
import dna_function.py



parser = argparse.ArgumentParser(description="Mix spectra data, sequencing data and Protein data base, to find new proteins")
parser.add_argument("spectra_file_name", type=str, help="the spectra file with format mzML")
parser.add_argument("sequencing_file_name", type=str, help="the fasta file, the best is to use Trinity with the fastq following by this software.\nThis soft is going to translate each sequence according to each reading frame.")
parser.add_argument("pdb_fasta_file", type=str, help="the fasta file which contains protein sequence")


a = parser.parse_args()






# the file with the mass data
spectra_file_name = a.spectra_file_name

# the file with the sequence after sequencing
fasta_seq_file = a.sequencing_file_name

# the fasta file which contains protein sequence
pdb_fasta_file = a.pdb_fasta_file




# We validate the msms with the sequencing data 
# We search the matchs between ms ms data and the fasta data
# we saved it in a csv file to use the data 

match_msms_vs_sequencing = threading.Thread(totalSearch, args=(spectra_file_name,fasta_seq_file,)

# we analyse the msms file with the pdb file         
match_msms_vs_PDB = threading.Thread(totalSearch, args=(spectra_file_name,pdb_fasta_file,)

# we analyse the alignement between the pdb file and the fasta file of the sequencing
# match_sequencing_vs_PDB = threading.Thread(, args=(fasta_seq_file,pdb_fasta_file,)


# we launch each thread process

threads_process = [match_msms_vs_sequencing, match_msms_vs_PDB] # , match_sequencing_vs_PDB 

for t in threads_process:
    t.start()

for t in threads_process:
    t.join()



# We analyse the previous results









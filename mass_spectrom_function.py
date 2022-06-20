from pyopenms import *
import threading

csv_file_spectra = threading.Semaphore()

list_of_tempory_files = [] 


def totalSearch(spectra_file_name, fasta_file_name):
    """
    To launch a search on a whole fasta file         

    """
    pass


def cut_a_fasta_file(fasta_file_name):
    """
    Cut a fasta file with several sequence to several fasta file with only one sequence


    """
    with open(fasta_file_name, 'r') as fasta_file:
        list_of_lines = fasta_file.readlines()
        # we delete the empty lines
        list_of_lines = [line for line in list_of_lines if line != "\n"] 
        index = 0 
        while index < len(list_of_lines)-1:
            # while it's not the end of file
            # we save the sequence and its head 
            head = list_of_lines[index]
            seq = list_of_lines[index+1]
            # we create a new temp file with this header and this sequence 
            with open("/tmp/"+fasta_file_name+str(index),'a') as short_fasta_file:
                short_fasta_file.write(head)
                short_fasta_file.write(seq)
            # we save the name of file to delete it at the end
            list_of_tempory_files.append("/tmp/"+fasta_file_name+str(index))
            index += 2
    return 0



def simpleSearch(spectra_file, fasta_file_name):
    """
    Peptide search example from https://pyopenms.readthedocs.io/en/latest/peptide_search.html

    Search the peptide from a spectra_file hitting with the sequence from a fasta_file_name




    """
    
    with open(fasta_file_name, 'r') as fasta_file:
        # the identification of the sequence
        main_key = fasta_file.readline()
        # the sequence
        sequence = fasta_file.readline()

    protein_ids = []
    peptide_ids = []

    SimpleSearchEngineAlgorithm().search(spectra_file, fasta_file, protein_ids, peptide_ids)
    

    pept = []
    
    # for each peptide hiting with the fasta sequence
    for peptide_id in peptide_ids:
        for hit in peptide_id.getHits():
            pept.append(hit.getSequence())
    

    csv_file_spectra.acquire()
    with open(fasta_file_name+"_vs_"+spectra_file+".csv",'a') as csv_file:
        for p in pept:
            csv_file.write(fasta_file_name+","+spectra_file+","+main_key+","+str(p)+"\n")
    csv_file_spectra.release()
    

    return protein_ids, peptide_ids




if __name__ == '__main__':
    import doctest
    # to test the different functions 
    doctest.testmod()
        
    cut_a_fasta_file("dna_function.py")





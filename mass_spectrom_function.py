from pyopenms import *
import threading

csv_file_spectra = threading.Semaphore()


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
    





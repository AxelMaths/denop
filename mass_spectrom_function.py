from pyopenms import *
import threading
import os


csv_file_spectra = threading.Semaphore()

list_of_tempory_files = [] 


def totalSearch(spectra_file_name, fasta_file_name):
    """
    To launch a search on a whole fasta file         

    """
    
    list_threadings=[]
     
    # we split the fasta file into multiple files
    cut_a_fasta_file(fasta_file_name)
    # we list the shorter files in tempory directory
    list_of_files_in_tmp = os.listdir(path="/tmp")
    # we select only the files about the current fasta_file_name
    list_of_files_for_the_current_fasta_file = [ "/tmp/"+a for a in list_of_files_in_tmp if fasta_file_name in a]
    # we display for the debug
    #printif(list_of_files_for_the_current_fasta_file)

    # for each file, we launch a thread to analyse each spectra file
    for short_fasta_file in list_of_files_for_the_current_fasta_file:
        # we create the thread, for the analysis
        list_threadings.append(threading.Thread(target=simpleSearch, args=(spectra_file, short_fasta_file,)))
        # we launch the thread
        list_threadings[-1].start()

    # we wait each job 
    for thread in list_threadings:
        thread.join()
     

    return 0



def cut_a_fasta_file(fasta_file_name):
    """
    Cut a fasta file with several sequence to several fasta file with only one sequence


    """
    with open(fasta_file_name, 'r') as fasta_file:
        fasta_file_name = fasta_file_name.split("/")[-1]
        print(fasta_file_name)
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
    
    # we launched the peptide tuto 
    protein_ids = []
    peptide_ids = []

    SimpleSearchEngineAlgorithm().search(spectra_file, fasta_file_name, protein_ids, peptide_ids)
    
    # we delete the path to use only the file name
    fasta_file_name = fasta_file_name.split("/")[-1]
    fasta_name = fasta_file_name.split(".")[0]
    print(fasta_name)

    pept = []
    
    # for each peptide hiting with the fasta sequence
    for peptide_id in peptide_ids:
        for hit in peptide_id.getHits():
            pept.append(hit.getSequence())
    
    # we saved the results in a csv file to be analysed
    csv_file_spectra.acquire()
    with open(fasta_file_name+"_vs_"+spectra_file+".csv",'a') as csv_file:
        for p in pept:
            csv_file.write(fasta_name+","+spectra_file+","+main_key+","+str(p)+"\n")
    csv_file_spectra.release()
    
    print(f"Done for the fasta file {fasta_file_name} and the spectra file {spectra_file}\n")

    return protein_ids, peptide_ids




if __name__ == '__main__':
    import doctest
    # to test the different functions 
    doctest.testmod()
        
    cut_a_fasta_file("dna_function.py")
    totalSearch("/Isiprod1/ext/Axel/Input/ms_data4/20022102.mgf","/Isiprod1/ext/Axel/Input/PDB/pdb_seqres.fasta")




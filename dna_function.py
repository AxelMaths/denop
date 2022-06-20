from Bio.Seq import Seq
import threading
from debug import printif

csv_amino_acid_semaphore = threading.Semaphore()


def seq_to_amin(seq, dico_to_save, key):
    """
    To translate a sequence of DNA to amino acid

    Warning : The sequence length has to be a multiple of 3 !


    >>> seq = "TTTATCGATCGATAA"
    >>> key = "key"
    >>> dico = dict()
    >>> seq_to_amin(seq, dico, key) 
    Seq('FIDR*')

    """
    dico_to_save
    dico_to_save[key] = Seq(seq).translate() 
    return dico_to_save[key]


def dna_to_amino_acid(sequence,titre):
    """
    Translate a sequence and its 6 reads frames.

    
    Complexity : O(len(sequence))
    
    To be with Threads
    
    >>> seq = "TTTATCTAA"
    >>> titre = "a title"
    >>> dna_to_amino_acid(seq,titre)
    {'a title+0': Seq('FI*'), 'a title-0': Seq('LDK'), 'a title+1': Seq('LS'), 'a title-1': Seq('*I'), 'a title+2': Seq('YL'), 'a title-2': Seq('R*')}
    >>> dna_to_amino_acid("",titre)
    {'a title+0': Seq(''), 'a title-0': Seq(''), 'a title+1': Seq(''), 'a title-1': Seq(''), 'a title+2': Seq(''), 'a title-2': Seq('')}

    """
    the_six_reads = dict()
    # With one reading sequence 
    
    # we set the sequence
    the_seq = Seq(sequence)#, generic_dna)
    printif(the_seq)
    # we set the reverse complementary sequence using Bio.Seq
    seq_comp = the_seq.reverse_complement()
    printif(seq_comp)
    
    length_seq = len(sequence)

    list_threads = []

    for i in [0,1,2]:
        # The following instructions use the threading module
        #
        # The length of the sequence are checked to be a multiple of 3, we cut the sequence to avoid exception.
        #
        # we prepare the translation of the i th frame for the starting sequence
        t1 = threading.Thread(target=seq_to_amin, args=(sequence[i:length_seq-(length_seq-i)%3], the_six_reads, titre+"+"+str(i), ))
        t1.start()
        list_threads.append(t1)

        # we prepare the translation of the i th frame for the reverse complementary sequence 
        t2 = threading.Thread(target=seq_to_amin, args=(str(seq_comp)[i:length_seq-(length_seq-i)%3], the_six_reads, titre+"-"+str(i), ))
        t2.start()
        list_threads.append(t2)

    for t in list_threads:
        t.join()
    
    printif(the_six_reads)
    return the_six_reads


def Fasta_to_amino_acid(fasta_file_name):
    """
    Create a csv file with the translated sequence from a fasta file

    The csv file is a database with the main key which is an integer.

    

    """
    list_threads = []
    with open(fasta_file_name, 'r', encoding = 'utf-8') as fIn:
        # we set the index for the database 
        index = 0 
        line  = fIn.readline()
        while(len(line) != 0):
            titre = line[:-1] 
            sequence = fIn.readline()[:-1]
           
            # we set the thread 
            thread_current = threading.Thread(target=write_fasta_to_amino_acid, args=(sequence, titre, fasta_file_name, index,) )
            thread_current.start()
            list_threads.append(thread_current)
            
            # we read the next line
            line = fIn.readline()
            # we increase the index of 6 
            index += 6
        
        # we wait each thread
        for t in list_threads:
            t.join()

    return index
            

def write_fasta_to_amino_acid(sequence, titre, fasta_file_name, index):
    """
    This function writes the six sequences of amino acid after the translation of the dna sequence.

    


    """
    translated_sequences_six_reading_frame = dna_to_amino_acid(sequence, titre)

    csv_amino_acid_semaphore.acquire() 
    with open(fasta_file_name+"_Translated_sequence.csv", 'a', encoding = 'utf-8') as fOut:
        for key in translated_sequences_six_reading_frame.keys():
            fOut.write(str(index)+" , "+str(key)+" , "+str(translated_sequences_six_reading_frame[key])+"\n")
            index+=1
    csv_amino_acid_semaphore.release()

    return index 





if __name__ == '__main__':
    import doctest
    # to test the different functions 
    doctest.testmod()

    # to test
    import os 
    os.system("touch fasta_test.fasta")
    with open("fasta_test.fasta",'a', encoding = 'utf-8') as f:
        f.write("> seq1\n")
        f.write("ATCGATCGATCGATCGATCG\n")
        f.write("> seq2\n")
        f.write("TCGATCAGTCAGCTATGCTAGCT\n")
    Fasta_to_amino_acid("fasta_test.fasta")
    
    os.system("rm fasta_test.fasta")


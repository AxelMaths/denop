from Bio.Seq import Seq
import threading

Debug = False

def printif(text, debug = Debug):
    """
    Print some data if and only if the variable debug is set True
    
    >>> printif("hello", debug=True)
    hello
    >>> printif("None", debug=False)
    
    """
    if debug:
        print(text)


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



if __name__ == '__main__':
    import doctest
    # to test the different functions 
    doctest.testmod()



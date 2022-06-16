from Bio.Seq import Seq
import threading

Debug = True

def printif(text, debug = Debug):
    """
    Print some data if and only if the variable debug is set True
    
    >>> printif("hello")
    hello
    >>> printif("None",debug=False)
    
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


if __name__ == '__main__':
    import doctest
    # to test the different functions 
    doctest.testmod()





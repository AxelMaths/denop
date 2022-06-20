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


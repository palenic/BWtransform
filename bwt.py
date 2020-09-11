# -*- coding: utf-8 -*-
""" Burrows-Wheeler transform

The Burrows-Wheeler transform is a very popular string manipulation 
algorithm at the basis of many bioinformatics tools as well as 
compression algorithms. This module provides functions for computing 
and inverting the Burrows-Wheeler transform on any input sequence 
given as a string.

Functions
--------
bwInverse(y)
    Invert the BW transform
bwTransform(x)
    Computes the BW transform
    
Notes
--------
Explanation of the algorithm used can be found in the help pages 
for the functions.

References
--------
.. [1] Burrows M, Wheeler D.J. A Block-Sorting Lossless 
Data Compression Algorithm (1994). Technical Report, 124, Digital 
Equipment Corporation.

Examples
--------
>>> bwInverse("g$actc")
'acctg'
>>> bwTransform("acctg")
'g$actc'
"""

import numpy as np
from re import findall 


class TerminatorError(ValueError):
    """
    This exception is raised if either bwTransform or bwInverse receive
    as input an invalid string.
    
    Attributes
    ----------
    msg : str
        Explanation of the error
    error_code: int
        An error code describing the error 
        
    """
    
    def __init__(self, msg, error_code):
        self.msg=msg
        self.error_code=error_code
        super().__init__(self.msg)
        

def bwTransform(x, terminator = "$"):
    """
    Compute the Burrows-Wheeler transform of an input string.

    Parameters
    ----------
    x : str
        The string to transform. Cannot contain the terminator character $.

    Raises
    ------
    TypeError
        If the input is not a string.
    TerminatorError
        If the input contains the terminator character ("$").

    Returns
    -------
    str
        The transformed string.
    
    Notes
    --------
    Implementation of the Burrows-Wheeler compression transformation, 
    originally described in [1]_, via construction of the suffix array, 
    as explained in section 2.2 of [2]_.
    The BWT of a string is obtained by definining a matrix of sorted 
    rotations of the string after adding a terminator character; 
    if the input is "hello":
    HELLO$       $HELLO
    ELLO$H       ELLO$H
    LLO$HE  -->  HELLO$
    LO$HEL       LLO$HE
    O$HELL       LO$HEL
    $HELLO       O$HELL
    The transform is the last column of the sorted matrix.
    The function adds a terminator $ to the string and derives a suffix 
    array S, which is an array containing the start position of all 
    suffixes sorted lexicographically. The BWT of X is defined as B[i]=$ 
    when S(i)=0 and B[i]=X[S(i)−1] otherwise; this is another way to 
    obtain the last column of the matrix without performing the rotations.
    
    References
    --------
    .. [1] Burrows M, Wheeler D.J. A Block-Sorting Lossless 
    Data Compression Algorithm (1994). Technical Report, 124, Digital 
    Equipment Corporation.
    .. [2] Li H, Durbin R. Fast and accurate short read alignment with 
    Burrows-Wheeler transform. Bioinformatics. 2009;25(14):1754-1760. 
    
    Examples
    --------
    >>> bwTransform("acctg")
    'g$actc'

    """
    ## input validation
    if not isinstance(x, str):
        raise TypeError("The input is not a string")
    if not isinstance(terminator, str):
        raise TerminatorError("The terminator must be a string of length 1",
                              error_code = 4)
    if len(terminator) != 1:
        raise TerminatorError("The terminator must be a string of length 1",
                              error_code = 5)
    if x.find(terminator) != -1:
        raise TerminatorError("The input cannot contain the terminator character ('{}'). Please provide another string.".format(terminator),
                              error_code = 1)
    ## adding the terminator
    allowed_ter = "$&*-%#"
    if allowed_ter.find(terminator) == -1:
        raise TerminatorError("Allowed characters for terminator are: $&*-%# ($ default)",
                              error_code = 6)
    x = x + terminator
    N = len(x)
    
    ## sorting all suffixes
    ## the list comprehension part creates the suffixes 
    ## and then np.argsort gives the indexes of the sorted
    ## suffixes
    suffix_array = np.argsort([x[i:] for i in range(N)])
    
    ## subtracting 1 from the index of sorted suffixes gives us the
    ## position on x of the last letter of the column. If 0,
    ## it becomes N-1 since the last character is $.
    bwt_indexes = np.where(suffix_array>0, suffix_array-1, N-1)
    
    ## building transformed string
    coded_seq = "".join([x[i] for i in bwt_indexes])
    return coded_seq

def bwInverse(y, terminator= "$"):
    """
    Invert the Burrows-Wheeler transform of a transformed input string.

    Parameters
    ----------
    y : str
        The transformed string. Must contain the terminator character $ 
        only once.

    Raises
    ------
    TypeError
        If the input is not a string.
    TerminatorError
        If the input does not contain the terminator character ("$") or 
        contains it more than once.

    Returns
    -------
    str
        The original string.
    
    Notes
    --------
    Implementation of the Burrows-Wheeler decompression as described 
    in [1]_, section 4.2. The transformed string y is the last column of
    the sorted matrix of rotations of the original string (L). We can 
    obtain the first column F by sorting lexicographically the characters 
    in L. Given how the matrix is built:
        1. L[i] cyclicly precedes F[i] in S for each i 
    To invert the transformation, we need to find a vector T that describes 
    the correspondence between each character in L and F such that: 
        2. F[T[j]]=L[j]
    Substituting i=T[j], from 1. we obtain:
        3. L[T[j]] precedes F[T[j]]
    From 2. we obtain:
        4. L[T[j]] cyclicly precedes L[j]
    Starting from the terminator character $=L[i], which is the last 
    character of S, L[T[i]] precedes it in S; then, substituting j=T[i], 
    L[T[j]] precedes T[j], and so on; this way, we reconstruct S back to 
    front.
    The function constructs an array P[N] and a dictionary C. C[ch] is the 
    total number of instances in L (the transformed string) of characters 
    preceding ch in the alphabet. P[i] is the number of instances of 
    character L[i] in the prefix L[0:i] of L. 
    T is defined as T[i]=P[i]+C[L[i]] for each i. In fact, the match of 
    L[i] in F can be found without building F by summing up
    a. how many characters in total precede L[i] in L (stored in C[L[i]]) 
    b. how many of the same character precede it in L[0:i] (stored in P[i])
    
    
    References
    --------
    .. [1] Burrows M, Wheeler D.J. A Block-Sorting Lossless 
    Data Compression Algorithm (1994). Technical Report, 124, Digital 
    Equipment Corporation.
    
    Examples
    --------
    >>> bwInverse("g$actc")
    'acctg'

    """

    ##input validation 
    if not isinstance(y, str):
        raise TypeError("The input is not a string")
    if not isinstance(terminator, str):
        raise TerminatorError("The terminator must be a string of length 1",
                              error_code = 4)
    if len(terminator) != 1:
        raise TerminatorError("The terminator must be a string of length 1",
                              error_code = 5)
    if y.find(terminator) == -1:
        raise TerminatorError("The input does not contain the terminator character ('{}'). Please provide a coded string.".format(terminator),
                              error_code = 2)
    if len(findall("\$", y)) > 1:
        raise TerminatorError("The terminator ('{}') cannot be present more than once in the coded string. Please provide another string.".format(terminator),
                              error_code = 3)   
    allowed_ter = "$&*-%#"
    if allowed_ter.find(terminator) == -1:
        raise TerminatorError("Allowed characters for terminator are: $&*-%# ($ default)",
                              error_code = 6)
    N = len(y)
    ## empty array for counts
    count_L = np.array(np.zeros(N))
    
    ## computing P array (which i call count_L)
    for i in range(N):
        char = y[i]
        char_occurrences = y[0:i].count(char)
        count_L[i] = char_occurrences
        
    ## computing C array (which i call count_alphabet and is a dictionary
    ## because I can use letters as keys)
    alphabet = np.sort(np.unique(list(y)))
    count_alphabet = {}
    ## initialize: no letter comes before the first one
    first_letter = alphabet[0]
    count_alphabet[first_letter] = 0
    
    ## could've done a dict comprehension but it was unreadable
    ## so I did a for loop
    for i in range(1, len(alphabet)):
        ## how many characters in L precede "letter" alphabetically?
        ## to do this I sum count_alphabet[prev_letter] and the counts 
        ## of "prev_letter" in L
        letter = alphabet[i]
        prev_letter = alphabet[i-1]
        count_alphabet[letter] = count_alphabet[prev_letter] + y.count(prev_letter)
    
    ## start from the terminator character and retrieve index of
    ## preceding character. Append them in a vector. Then find character
    ## preceding that one and so on. This gives us the letters of S in 
    ## inverse order (as their position on L)
    i = y.find(terminator)
    decoded_indexes = np.array(np.zeros(N), dtype="int")
    
    for j in range(N-1, -1, -1):
        decoded_indexes[j] = i
        i = int(count_L[i] + count_alphabet[y[i]])
        
    ## building original string + omitting terminator
    decoded_seq="".join([y[i] for i in decoded_indexes])[0:N-1]    
    ## return the string 
    return decoded_seq

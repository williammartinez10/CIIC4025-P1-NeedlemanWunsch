
''' 
Implementation of Needleman-Wunsch algorithm, using Dynamic Programming

    This algorithm tries to align an entire sequence, align all letters
    from query and target, and it is suitable for closely related sequences.

    This tool accepts comma separated value files (.csv).

    To be able to use it, this script requires 'numpy' to be installed within the Python
    environment you are running it in.

@author: William A. Martínez Martínez

'''

import numpy as np
import sys


def S(A, B):
    '''
        Scoring Matrix

        Parameters
        ----------
        A : str
            Character in first biological sequence
        B : str
            Character in second biological sequence

        Returns
        -------
        int
            Match reward if equal, mismatch penalty if not equal, gap penalty otherwise

    '''
    if A == B:
        return matchReward
    elif A != B:
        return mismatchPenalty
    else:
        return gapPenalty


def needleman_wunsch(S1, S2):
    '''
        Algorithm for global alignment of two biological sequences.

        Parameters
        ----------
        S1 : str
            First biological sequence
        S2 : str
            Second biological sequence

        Returns
        -------
        tuple
            Aligned biological sequences 1 and 2, and alignment score
    '''
    m = len(S1)
    n = len(S2)
    F = np.zeros(( m+1 , n+1 )).astype(int) # Creates an m x n matrix for two given sequences
    
    # 1 - Initialization Step : Initializes matrix with posible scores
    for j in range(n+1):
        F[0][j] = gapPenalty * j
    for i in range(m+1):
        F[i][0] = gapPenalty * i

    # 2 - Matrix Fill Step : Fills matrix with maximum scores
    for i in range(1, m+1):
        for j in range(1, n+1):
            F[i][j] = max(  F[i-1][j-1] + S(S1[i-1], S2[j-1]), 
                            F[i][j-1] + gapPenalty, 
                            F[i-1][j] + gapPenalty
                        )

    # 3 - Backtracking Step : Traceback for appropriate alignment
    alignedS1 = ""
    alignedS2 = ""
    while m > 0 or n > 0:
        if F[m][n] == ( F[m-1][n] + gapPenalty ): # Checks left element of current position in matrix for comparing
            alignedS1 = S1[m-1] + alignedS1
            alignedS2 = "-" + alignedS2
            m -= 1
        elif F[m][n] == ( F[m][n-1] + gapPenalty ): # Checks upper element of current position in matrix for comparing
            alignedS1 = "-" + alignedS1
            alignedS2 = S2[n-1] + alignedS2
            n -= 1
        else: # Accepts diagonal element of current position in matrix
            alignedS1 = S1[m-1] + alignedS1
            alignedS2 = S2[n-1] + alignedS2
            m -= 1
            n -= 1

    alignmentScore = F[-1][-1]
    
    return (alignedS1, alignedS2, alignmentScore)


###
#   Main Program
###
if len(sys.argv) > 1:
    if __name__ == '__main__':
        # Reward and penalty system
        matchReward = 1
        mismatchPenalty = -1
        gapPenalty = -2

        # In command line, run "python3 ./<PythonFilename>.py <CSVFilename>.csv"
        # Parses input file
        with open(sys.argv[1], mode='r', encoding='utf-8-sig') as f:
            for line in f:
                fields = line.strip().split(',') # Separates rows and removes newline character from the last element in them
                sequence1 = fields[0]
                sequence2 = fields[1]
                if not any(map(str.isdigit, sequence1 or sequence2)): # Checks if a sequence is valid
                    result = needleman_wunsch(sequence1, sequence2)
                    print(' '.join(map(str, result)))
                else:
                    continue

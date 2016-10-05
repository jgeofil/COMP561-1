#!/usr/bin/python

import sys, getopt
import numpy as np

def main(argv):
    helpMsg = 'test.py -i <inputfile> -a <multipleOfThreeScore> -b <nonMultipleOfThreeScore> -m <scoreForMatches> -t <transitionScore> -v <transversionScore>'
    try:
        opts, args = getopt.getopt(argv,"hi:a:b:m:t:v:")
    except getopt.GetoptError:
        print(helpMsg)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(helpMsg)
            sys.exit()
        elif opt == "-i":
            inputfile = arg
        elif opt == "-a":
            GAP_PEN_A = int(arg)
        elif opt == "-b":
            GAP_PEN_B = int(arg)
        elif opt == "-m":
            MATCH_SCORE = int(arg)
        elif opt == "-t":
            TRANSITION_SCORE = int(arg)
        elif opt == "-v":
            TRANSVERSION_SCORE = int(arg)

    with open(inputfile) as f:
        count = 0
        SEQ1 = ''
        SEQ2 = ''
        for line in f:
            if line.startswith('>'):
                count += 1
            elif count == 1:
                SEQ1 += line
            elif count == 2:
                SEQ2 += line
        f.close()
        SEQ1 = SEQ1.rstrip()
        SEQ2 = SEQ2.rstrip()

    M = len(SEQ1)
    M1 = M + 1
    N = len(SEQ2)
    N1 = N + 1

    D = np.empty([M1,N1])
    VA1, VA2, VA3 = np.empty([M1,N1]), np.empty([M1,N1]), np.empty([M1,N1])
    HA1, HA2, HA3 = np.empty([M1,N1]), np.empty([M1,N1]), np.empty([M1,N1])

    PUR = ('A', 'G')
    PYR = ('T', 'C')
    BASES = PUR + PYR

    SCORING_M = {x: {y: MATCH_SCORE if x==y else TRANSITION_SCORE if (x in PUR and y in PUR) or (x in PYR and y in PYR) else TRANSVERSION_SCORE for y in BASES} for x in BASES}

    # Initialisation
    for i in range(0,M1):
        D[i,0] = i * GAP_PEN_B
        HA1[i,0], HA2[i,0], HA3[i,0] = float("-inf"), float("-inf"), float("-inf")
        VA1[i,0], VA2[i,0], VA3[i,0] = float("-inf"), float("-inf"), float("-inf")
        VA3[i,0] = GAP_PEN_A*i if i%3 == 0 else float("-inf")
        VA1[i,0] = GAP_PEN_A*i if i%1 == 0 else float("-inf")
        VA2[i,0] = GAP_PEN_A*i if i%2 == 0 else float("-inf")

    for j in range(0,N1):
        D[0,j] = j * GAP_PEN_B
        VA1[0,j], VA2[0,j], VA3[0,j] = float("-inf"), float("-inf"), float("-inf")
        HA1[0,j], HA2[0,j], HA3[0,j] = float("-inf"), float("-inf"), float("-inf")
        HA3[0,j] = GAP_PEN_A*j if j%3 == 0 else float("-inf")
        HA1[0,j] = GAP_PEN_A*j if j%1 == 0 else float("-inf")
        HA2[0,j] = GAP_PEN_A*j if j%2 == 0 else float("-inf")

    print('Recursion...')
    # Recursion
    for i in range(1,M1):
        for j in range(1,N1):

            match = SCORING_M[SEQ1[i-1]][SEQ2[j-1]]
            dij1 = D[i, j-1]
            di1j = D[i-1, j]

            Ddiag = D[i-1, j-1] + match
            DtopB = di1j + GAP_PEN_B
            DleftB = dij1 + GAP_PEN_B
            Vtop = VA3[i-1,j-1] + match
            Hleft = HA3[i-1,j-1] + match

            D[i,j] = max(DtopB, DleftB, Vtop, Hleft, Ddiag)

            DtopA = di1j + GAP_PEN_A
            VA3topA = VA3[i-1, j] + GAP_PEN_A

            VA1[i,j] = max(DtopA, VA3topA)
            VA2[i,j] = VA1[i-1,j] + GAP_PEN_A
            VA3[i,j] = VA2[i-1,j] + GAP_PEN_A

            DleftA = dij1 + GAP_PEN_A
            HA3leftA = HA3[i, j-1] + GAP_PEN_A

            HA1[i,j] = max(DleftA, HA3leftA)
            HA2[i,j] = HA1[i, j-1] + GAP_PEN_A
            HA3[i,j] = HA2[i, j-1] + GAP_PEN_A
    #HA1[HA1>1000] = 100

    for row in HA1[7:11, 68:72]:
        row[row<-10000] =10
    for row in HA2[7:11, 68:72]:
        row[row<-10000] =10
    for row in HA3[7:11, 68:72]:
        row[row<-10000] =10
    print(D.astype('int')[7:11, 68:72])
    print(HA1.astype('int')[7:11, 68:72])
    print(HA2.astype('int')[7:11, 68:72])
    print(HA3.astype('int')[7:11, 68:72])

    # Traceback
    print('Tracing back...')
    ali1, ali2 = '',''
    i, j = M, N
    state = 'D' if D[i,j] > HA3[i,j] and D[i,j] > VA3[i,j] else 'H3' if HA3[i,j] > VA3[i,j] else 'V3'
    score = D[i,j] if state == 'D' else HA3[i,j] if state == 'H3' else VA3[i,j]
    print(score)

    while i > 0 or j > 0:
        print(state)
        print(i)
        print(j)

        if state == 'D':


            match = SCORING_M[SEQ1[i-1]][SEQ2[j-1]]

            scoreDiag = D[i-1, j-1] + match if i > 0 and j > 0 else float('-inf')
            scoreTop = D[i-1, j] + GAP_PEN_B
            scoreLeft = D[i, j-1] + GAP_PEN_B
            scoreH3 = HA3[i-1, j-1] + match
            scoreV3 = VA3[i-1, j-1] + match

            print(scoreDiag)
            print(scoreTop)
            print(scoreLeft)
            print(scoreH3)
            print(scoreV3)
            print(score)



            if score == scoreTop:
                ali1 += SEQ1[i-1]
                ali2 += '-'
                state = 'D'
                i -= 1
                score = D[i,j]

            elif score == scoreLeft:
                ali1 += '-'
                ali2 += SEQ2[j-1]
                state = 'D'
                j -= 1
                score = D[i,j]

            elif score == scoreH3:
                ali1 += SEQ1[i-1]
                ali2 += SEQ2[j-1]
                state = 'H3'
                i -= 1
                j -= 1
                score = HA3[i,j]

            elif score == scoreV3:
                ali1 += SEQ1[i-1]
                ali2 += SEQ2[j-1]
                state = 'V3'
                i -= 1
                j -= 1
                score = VA3[i,j]
            elif score == scoreDiag:
                ali1 += SEQ1[i-1]
                ali2 += SEQ2[j-1]
                state = 'D'
                i -= 1
                j -= 1
                score = D[i,j]

        elif state == 'H1':
            if i == 0:
                ali1 += '-'
                ali2 += SEQ2[j-1]
                state = 'H3'
                j -= 1
                score = HA3[i,j]
            else:
                scoreH3 = HA3[i, j-1]
                scoreLeft = D[i, j-1]

                if score == scoreH3 + GAP_PEN_A:
                    ali1 += '-'
                    ali2 += SEQ2[j-1]
                    state = 'H3'
                    j -= 1
                    score = HA3[i,j]
                elif score == scoreLeft + GAP_PEN_A:
                    ali1 += '-'
                    ali2 += SEQ2[j-1]
                    state = 'D'
                    j -= 1
                    score = D[i,j]
        elif state == 'H2':
            ali1 += '-'
            ali2 += SEQ2[j-1]
            state = 'H1'
            j -= 1
            score = HA1[i,j]
        elif state == 'H3':
            ali1 += '-'
            ali2 += SEQ2[j-1]
            state = 'H2'
            j -= 1
            score = HA2[i,j]
        elif state == 'V1':
            if j == 0:
                ali1 += SEQ1[i-1]
                ali2 += '-'
                state = 'V3'
                i -= 1
                score = VA3[i,j]
            else:
                scoreV3 = VA3[i-1, j]
                scoreTop = D[i-1, j]

                if score == scoreV3 + GAP_PEN_A:
                    ali1 += SEQ1[i-1]
                    ali2 += '-'
                    state = 'V3'
                    i -= 1
                    score = VA3[i,j]
                elif score == scoreTop + GAP_PEN_A:
                    ali1 += SEQ1[i-1]
                    ali2 += '-'
                    state = 'D'
                    i -= 1
                    score = D[i,j]
        elif state == 'V2':
            ali1 += SEQ1[i-1]
            ali2 += '-'
            state = 'V1'
            i -= 1
            score = VA1[i,j]
        elif state == 'V3':
            ali1 += SEQ1[i-1]
            ali2 += '-'
            state = 'V2'
            i -= 1
            score = VA2[i,j]

    def printAlignement(a,b):
        n = 80
        a = [a[i:i+n] for i in range(0, len(a), n)]
        b = [b[i:i+n] for i in range(0, len(b), n)]

        for i in range(0, len(a)):
            head = '|    '*(n/5)
            print(i*80+1)
            print(head)
            print(a[i])
            print(b[i])
            print(' ')


    printAlignement(ali1[::-1],ali2[::-1])

if __name__ == "__main__":
   main(sys.argv[1:])

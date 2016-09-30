import numpy as np

SEQ1 = 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTGCAGGAAACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACTGTGAGAACTCTGAGGACAAAGCAGCGGATACAACCTCAAAAGACGTCTGTCTACATTGAATTGGGATCTGATTCTTCTGAAGATACCGTTAATAAGGCAACTTATTGCAGTGTGGGAGATCAAGAATTGTTACAAATCACCCCTCAAGGAACCAGGGATGAAATCAGTTTGGATTCTGCAAAAAAGGCTGCTTGTGAATTTTCTGAGACGGATGTAACAAATACTGAACATCATCAACCCAGTAATAATGATTTGAACACCACTGAGAAGCGTGCAGCTGAGAGGCATCCAGAAAAGTATCAGGGTAGTTCTGTTTCAAACTTGCATGTGGAGCCATGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACTAAAGACAGAATGAATGTAGAAAAGGCTGAATTCTGTAATAAAAGCAAACAGCCTGGCTTAGCAAGGAGCCAACATAACAGATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACTCCCAGCACAGAAAAAAAGGTAGATCTGAATGCTGATCCCCTGTGTGAGAGAAAAGAATGGAATAAGCAGAAACTGCCATGCTCAGAGAATCCTAGAGATACTGAAGATGTTCCTTGGATAACACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATTGGACGTTCTAAATGAGGTAGATGAATATTCTGGTTCTTCAGAGAAAATAGACTTACTGGCCAGTGATCCTCATGAGGCTTTAATATGTAAAAGTGAAAGAGTTCACTCCAAATCAGTAGAGAGTAATATTGAAGACAAAATATTTGGGAAAACCTATCGGAAGAAGGCAAGCCTCCCCAACTTAAGCCATGTAACTGAAAATCTAATTATAGGAGCATTTGTTACTGAGCCACAGATAATACAAGAGCGTCCCCTCACAAATAAATTAAAGCGTAAAAGGAGACCTACATCAGGCCTTCATCCTGAGGATTTTATCAAGAAAGCAGATTTGGCAGTTCAAAAGACTCCTGAAATGATAAATCAGGGAACTAACCAAACGGAGCAGAATGGTCAAGTGATGAATATTACTAATAGTGGTCATGAGAATAAAACAAAAGGTGATTCTATTCAGAATGAGAAAAATCCTAACCCAATAGAATCACTCGAAAAAGAATCTGCTTTCAAAACGAAAGCTGAACCTATAAGCAGCAGTATAAGCAATATGGAACTCGAATTAAATATCCACAATTCAAAAGCACCTAAAAAGAATAGGCTGAGGAGGAAGTCTTCTACCAGGCATATTCATGCGCTTGAACTAGTAGTCAGTAGAAATCTAAGCCCACCTAATTGTACTGAATTGCAAATTGATAGTTGTTCTAGCAGTGAAGAGATAAAGAAAAAAAAGTACAACCAAATGCCAGTCAGGCACAGCAGAAACCTACAACTCATGGAAGGTAAAGAACCTGCAACTGGAGCCAAGAAGAGTAACAAGCCAAATGAACAGACAAGTAAAAGACATGACAGCGATACTTTCCCAGAGCTGAAGTTAACAAATGCACCTGGTTCTTTTACTAAGTGTTCAAATACCAGTGAACTTAAAGAATTTGTCAATCCTAGCCTTCCAAGAGAAGAAAAAGAAGAGAAACTAGAAACAGTTAAAGTGTCTAATAATGCTGAAGACCCCAAAGATCTCATGTTAAGTGGAGAAAGGGTTTTGCAAACTGAAAGATCTGTAGAGAGTAGCAGTATTTCATTGGTACCTGGTACTGATTATGGCACTCAGGAAAGTATCTCGTTACTGGAAGTTAGCACTCTAGGGAAGGCAAAAACAGAACCAAATAAATGTGTGAGTCAGTGTGCAGCATTTGAAAACCCCAAGGGACTAATTCATGGTTGTTCCAAAGATAATAGAAATGACACAGAAGGCTTTAAGTATCCATTGGGACATGAAGTTAACCACAGTCGGGAAACAAGCATAGAAATGGAAGAAAGTGAACTTGATGCTCAGTATTTGCAGAATACATTCAAGGTTTCAAAGCGCCAGTCATTTGCTCCGTTTTCAAATCCAGGAAATGCAGAAGAGGAATGTGCAACATTCTCTGCCCACTCTGGGTCCTTAAAGAAACAAAGTCCAAAAGTCACTTTTGAATGTGAACAAAAGGAAGAAAATCAAGGAAAGAATGAGTCTAATATCAAGCCTGTACAGACAGTTAATATCACTGCAGGCTTTCCTGTGGTTGGTCAGAAAGATAAGCCAGTTGATAATGCCAAATGTAGTATCAAAGGAGGCTCTAGGTTTTGTCTATCATCTCAGTTCAGAGGCAACGAAACTGGACTCATTACTCCAAATAAACATGGACTTTTACAAAACCCATATCGTATACCACCACTTTTTCCCATCAAGTCATTTGTTAAAACTAAATGTAAGAAAAATCTGCTAGAGGAAAACTTTGAGGAACATTCAATGTCACCTGAAAGAGAAATGGGAAATGAGAACATTCCAAGTACAGTGAGCACAATTAGCCGTAATAACATTAGAGAAAATGTTTTTAAAGAAGCCAGCTCAAGCAATATTAATGAAGTAGGTTCCAGTACTAATGAAGTGGGCTCCAGTATTAATGAAATAGGTTCCAGTGATGAAAACATTCAAGCAGAACTAGGTAGAAACAGAGGGCCAAAATTGAATGCTATGCTTAGATTAGGGGTTTTGCAACCTGAGGTCTATAAACAAAGTCTTCCTGGAAGTAATTGTAAGCATCCTGAAATAAAAAAGCAAGAATATGAAGAAGTAGTTCAGACTGTTAATACAGATTTCTCTCCATATCTGATTTCAGATAACTTAGAACAGCCTATGGGAAGTAGTCATGCATCTCAGGTTTGTTCTGAGACACCTGATGACCTGTTAGATGATGGTGAAATAAAGGAAGATACTAGTTTTGCTGAAAATGACATTAAGGAAAGTTCTGCTGTTTTTAGCAAAAGCGTCCAGAAAGGAGAGCTTAGCAGGAGTCCTAGCCCTTTCACCCATACACATTTGGCTCAGGGTTACCGAAGAGGGGCCAAGAAATTAGAGTCCTCAGAAGAGAACTTATCTAGTGAGGATGAAGAGCTTCCCTGCTTCCAACACTTGTTATTTGGTAAAGTAAACAATATACCTTCTCAGTCTACTAGGCATAGCACCGTTGCTACCGAGTGTCTGTCTAAGAACACAGAGGAGAATTTATTATCATTGAAGAATAGCTTAAATGACTGCAGTAACCAGGTAATATTGGCAAAGGCATCTCAGGAACATCACCTTAGTGAGGAAACAAAATGTTCTGCTAGCTTGTTTTCTTCACAGTGCAGTGAATTGGAAGACTTGACTGCAAATACAAACACCCAGGATCCTTTCTTGATTGGTTCTTCCAAACAAATGAGGCATCAGTCTGAAAGCCAGGGAGTTGGTCTGAGTGACAAGGAATTGGTTTCAGATGATGAAGAAAGAGGAACGGGCTTGGAAGAAAATAATCAAGAAGAGCAAAGCATGGATTCAAACTTAGGTGAAGCAGCATCTGGGTGTGAGAGTGAAACAAGCGTCTCTGAAGACTGCTCAGGGCTATCCTCTCAGAGTGACATTTTAACCACTCAGCAGAGGGATACCATGCAACATAACCTGATAAAGCTCCAGCAGGAAATGGCTGAACTAGAAGCTGTGTTAGAACAGCATGGGAGCCAGCCTTCTAACAGCTACCCTTCCATCATAAGTGACTCTTCTGCCCTTGAGGACCTGCGAAATCCAGAACAAAGCACATCAGAAAAAGCAGTATTAACTTCACAGAAAAGTAGTGAATACCCTATAAGCCAGAATCCAGAAGGCCTTTCTGCTGACAAGTTTGAGGTGTCTGCAGATAGTTCTACCAGTAAAAATAAAGAACCAGGAGTGGAAAGGTCATCCCCTTCTAAATGCCCATCATTAGATGATAGGTGGTACATGCACAGTTGCTCTGGGAGTCTTCAGAATAGAAACTACCCATCTCAAGAGGAGCTCATTAAGGTTGTTGATGTGGAGGAGCAACAGCTGGAAGAGTCTGGGCCACACGATTTGACGGAAACATCTTACTTGCCAAGGCAAGATCTAGAGGGAACCCCTTACCTGGAATCTGGAATCAGCCTCTTCTCTGATGACCCTGAATCTGATCCTTCTGAAGACAGAGCCCCAGAGTCAGCTCGTGTTGGCAACATACCATCTTCAACCTCTGCATTGAAAGTTCCCCAATTGAAAGTTGCAGAATCTGCCCAGAGTCCAGCTGCTGCTCATACTACTGATACTGCTGGGTATAATGCAATGGAAGAAAGTGTGAGCAGGGAGAAGCCAGAATTGACAGCTTCAACAGAAAGGGTCAACAAAAGAATGTCCATGGTGGTGTCTGGCCTGACCCCAGAAGAATTTATGCTCGTGTACAAGTTTGCCAGAAAACACCACATCACTTTAACTAATCTAATTACTGAAGAGACTACTCATGTTGTTATGAAAACAGATGCTGAGTTTGTGTGTGAACGGACACTGAAATATTTTCTAGGAATTGCGGGAGGAAAATGGGTAGTTAGCTATTTCTGGGTGACCCAGTCTATTAAAGAAAGAAAAATGCTGAATGAGCATGATTTTGAAGTCAGAGGAGATGTGGTCAATGGAAGAAACCACCAAGGTCCAAAGCGAGCAAGAGAATCCCAGGACAGAAAGATCTTCAGGGGGCTAGAAATCTGTTGCTATGGGCCCTTCACCAACATGCCCACAGATCAACTGGAATGGATGGTACAGCTGTGTGGTGCTTCTGTGGTGAAGGAGCTTTCATCATTCACCCTTGGCACAGGTGTCCACCCAATTGTGGTTGTGCAGCCAGATGCCTGGACAGAGGACAATGGCTTCCATGCAATTGGGCAGATGTGTGAGGCACCTGTGGTGACCCGAGAGTGGGTGTTGGACAGTGTAGCACTCTACCAGTGCCAGGAGCTGGACACCTACCTGATACCCCAGATCCCCCACAGCCACTACTGA'
SEQ2 = 'ATGGATTTATCTGCCGTCCAAATTCAAGAAGTACAAAATGTCCTTCATGCTATGCAGAAAATCTTAGAGTGTCCGATCTGTTTGGAACTGATCAAAGAACCTGTTTCCACAAAGTGTGACCACATATTTTGCAAATTTTGTATGCTGAAACTTCTTAACCAGAAGAAAGGGCCTTCACAATGTCCTTTGTGTAAGAATGAGATAACCAAAAGGAGCCTACAGGGAAGCACAAGGTTTAGTCAGCTTGCTGAAGAGCTGCTGAGAATAATGGCTGCTTTTGAGCTTGACACGGGAATGCAGCTTACAAATGGTTTTAGTTTTTCAAAAAAGAGAAATAATTCTTGTGAGCGTTTGAATGAGGAGGCGTCGATCATCCAGAGCGTGGGCTACCGGAACCGTGTCAGAAGGCTTCCCCAGGTCGAACCTGGAAATGCCACCTTGAAGGACAGCCTAGGTGTCCAGCTGTCTAACCTTGGAATCGTGAGATCAGTGAAGAAAAACAGGCAGACCCAACCTCGAAAGAAATCTGTCTACATTGAACTAGACTCTGATTCTTCTGAAGAGACAGTAACTAAGCCAGGTGATTGCAGTGTGAGAGACCAGGAATTGTTACAGACCGCCCCTCAAGAAGCTGGAGATGAAGGCAAGCTGCACTCTGCAGAAGAGGCTGCTTGTGAGTTTTCTGAGGGCATAAGAAACATTGAACATCATCAATGCAGTGATGATTTAAACCCTACTGAGAATCATGCAACTGAAAGGCATCCAGAAAAATGTCAGAGTATTTCTATTTCAAATGTGTGTGTGGAGCCATGTGGCACAGATGCTCATGCCAGCTCATTACAGCCTGAGACCAGCAGTTTATTGCTCATTGAAGACAGAATGAATGCAGAAAAGGCTGAATTCTGTAATAAAAGCAAACAGCCTGGCATAGCAGTGAGCCAGCAGAGCAGATGGGCTGCAAGTAAAGGAACATGTAACGACAGGCAGGTTCCCAGCACTGGGGAAAAGGTAGGTCCAAACGCTGACTCCCTTAGTGATAGAGAGAAGTGGACTCACCCGCAAAGTCTGTGCCCTGAGAATTCTGGAGCTACCACCGATGTTCCTTGGATAACACTAAATAGCAGCGTTCAGAAAGTTAATGAGTGGTTTTCCAGAACTGGTGAAATGTTAACTTCTGACAGCGCATCTGCCAGGAGGCACGAGTCAAATGCTGAAGCAGCTGTTGTGTTGGAAGTTTCAAACGAAGTGGATGGGGGTTTTAGTTCTTCAAGGAAAACAGACTTAGTAACCCCCGACCCCCATCATACTTTAATGTGTAAAAGTGGAAGAGACTTCTCCAAACCAGTAGAGGATAATATCAGTGATAAAATATTTGGGAAATCCTATCAGAGAAAGGGAAGCCGCCCTCACCTGAACCATGTGACTGAAATTATAGGCACATTTATTACAGAACCACAGATAACACAAGAGCAGCCCTTCACAAATAAATTAAAACGTAAGAGAAGTACATCCCTTCAACCTGAGGACTTCATCAAGAAAGCAGATTCAGCAGGTGTTCAAAGGACTCCTGACAACATAAATCAGGGAACTGACCTAATGGAGCCAAATGAGCAAGCAGTGAGTACTACCAGTAACTGTCAGGAGAACAAAATAGCAGGTAGTAATCTCCAGAAAGAGAAAAGCGCTCATCCAACTGAATCATTGAGAAAGGAACCTGCTTCCACAGCAGGAGCCAAATCTATAAGCAACAGTGTAAGTGATTTGGAGGTAGAATTAAACGTCCACAGTTCAAAAGCACCTAAGAAAAATAGGCTGAGGAGGAAGTCTTCTATCAGGTGTGCTCTTCCACTTGAACCAATCAGTAGAAATCCAAGCCCACCTACTTGTGCTGAGCTTCAAATCGATAGTTGTGGTAGCAGTGAAGAAACAAAGAAAAACCATTCCAACCAACAGCCAGCCGGGCACCTTAGAGAGCCTCAACTCATCGAAGACACTGAACCTGCAGCGGATGCCAAGAAGAACGAGCCAAATGAACACATAAGGAAGAGACGTGCCAGCGATGCTTTCCCAGAAGAGAAATTAATGAACAAAGCTGGTTTATTAACTAGCTGTTCAAGTCCTAGAAAATCTCAAGGGCCTGTCAATCCCAGCCCTCAGAGAACAGGAACAGAGCAACTTGAAACACGCCAAATGTCTGACAGTGCCAAAGAACTCGGGGATCGGGTCCTAGGAGGAGAGCCCAGTGGCAAAACCACTGACCGATCTGAGGAGAGCACCAGCGTATCCTTGGTATCTGACACTGACTACGACACTCAGAACAGTGTCTCAGTCCTGGACGCTCACACTGTCAGATATGCAAGAACAGGATCCGCTCAGTGTATGACTCAGTTTGTAGCAAGCGAAAACCCCAAGGAACTCGTCCATGGCTCTAACAATGCTGGGAGTGGCACAGAGGGTCTCAAGCCCCCCTTGAGACACGCGCTTAACCTCAGTCAGGAGAAAGTAGAAATGGAAGACAGTGAACTTGATACTCAGTATTTGCAGAATACATTTCAAGTTTCAAAGCGTCAGTCATTTGCTTTATTTTCAAAACCTAGAAGTCCCCAAAAGGACTGTGCTCACTCTGTGCCCTCAAAGGAACTGAGTCCAAAGGTGACAGCTAAAGGTAAACAAAAAGAACGTCAGGGACAGGAAGAATTTGAAATCAGTCACGTACAAGCAGTTGCGGCCACAGTGGGCTTACCTGTGCCCTGTCAAGAAGGTAAGCTAGCTGCTGATACAATGTGTGATAGAGGTTGTAGGCTTTGTCCATCATCTCATTACAGAAGCGGGGAGAATGGACTCAGCGCCACAGGTAAATCAGGAATTTCACAAAACTCACATTTTAAACAATCAGTTTCTCCCATCAGGTCATCTATAAAAACTGACAATAGGAAACCTCTGACAGAGGGACGATTTGAGAGACATACATCATCAACTGAGATGGCGGTGGGAAATGAGAACATTCTTCAGAGTACAGTGCACACAGTTAGCCTGAATAACAGAGGAAATGCTTGTCAAGAAGCCGGCTCGGGCAGTATTCATGAAGTATGTTCCACTGGTGACTCCTTCCCAGGACAACTAGGTAGAAACAGAGGGCCTAAGGTGAACACTGTGCCTCCATTAGATAGTATGCAGCCTGGTGTCTGTCAGCAAAGTGTTCCTGTAAGTGATAAGTATCTTGAAATAAAAAAGCAGGAGGGTGAGGCTGTCTGTGCAGACTTCTCTCCATGTCTATTCTCAGACCATCTTGAGCAATCTATGAGTGGTAAGGTTTTTCAGGTTTGCTCTGAGACACCTGATGACCTGCTGGATGATGTTGAAATACAGGGACATACTAGCTTTGGTGAAGGTGACATAATGGAGAGATCTGCTGTCTTTAACGGAAGCATCCTGAGAAGGGAGTCCAGTAGGAGCCCTAGTCCTGTAACCCATGCATCGAAGTCTCAGAGTCTCCACAGAGCGTCTAGGAAATTAGAATCGTCAGAAGAGAGCGACTCCACTGAGGATGAAGATCTTCCCTGCTTCCAACACTTACTGAGCAGAATAAGCAACACACCTGAGCTTACCAGATGCAGCAGTGCTGTGACACAGCGTATGCCAGAGAAAGCGGAGGGGACCCAAGCACCATGGAAGGGTAGCAGCAGTGACTGCAATAATGAGGTGATCATGATAGAGGCATCTCAGGAGCATCAGTTTAGTGAGGATCCAAGATGCTCTGGCAGCATGTTCTCTTCACAGCACAGTGCTGCCCAAGGGTCAACTGCAAATGCAAACTCCCAGGATTCCAATTTTATTCCACCTTCCAAACAGAGGAGTCACCAGTGTGGGAATGAGGAAGCTTTCCTAAGTGACAAGGAATTGATTTCAGATAACGAGGAAATGGCAACTTGCCTAGAAGAGGATAATGACCAAGAAGAGGATAGTATAATCCCAGATTCAGAGGCATCCGGATACGAGAGTGAAACAAACCTTTCTGAAGACTGCTCGCAGAGTGATATTTTAACCACTCAGCAGCGGGCGACCATGAAGTATAACCTGATAAAGCTGCAGCAGGAAATGGCTCACCTGGAAGCTGTGCTGGAGCAGCGTGGGAACCAGCCTTCTGGCCACTCCCCTTCCCTCCTAGCGGACCCTTGTGCCCTGGAAGACCTGCCAGATCTGGAACCAAACATGTCAGGAGCAGCAATTTTAACTTCAAAGAACATTAATGAGAATCCTGTAAGCCAAAATTTGAAGAGCGCTTGTGATGACAAATTCCAACTACAACATCTGGAGGGTCCCACCAGTGGAGATGACGAGTCAGGGATGGGAAGGCCTTCCCCTTTTAAATCTCCGTTGGCAGGCAGTAGGGGCTCTGCACATGGCTGCTCTAGGCATCTTCAAAAGAGAAACTCCCCCTCTCAGGAGGAGCTCCTCCAGCCTGCTGGATCAGAGGCGTCATCTGAGCCACACAATTCAACAGGGCAGTCTTGCCTGCCAAGGCGAGAGCTAGAAGGAACCCCATACCTGGGATCTGGAATCAGCCTTTTCTCTAGTAGAGACCCCGAATCTGAGTCCCCTAAAGAGCCAGCCCACATTGGCACCACACCAGCTTCAACCTCTGCACTGAAAATACCCCAAGGTCAAGTTGCTTTCCGGAGTGCAGCTGCTGCTGGTGCTGATAAAGCAGTGGTAGGAATTGTGAGCAAGATAAAGCCGGAATTGACATCTTCAGAAGAAAGAGCGGATAGAGACATATCCATGGTGGTGTCAGGCTTGACCCCCAAAGAAGTAATGACCGTGCAAAAGTTTGCTGAAAAATACCGCCTCACTTTAACTGACGCAATTACTGAGGAGACTACACATGTAATTATAAAAACAGATGCGGAGTTTGTGTGTGAGCGGACACTGAAATATTTTCTGGGCATTGCAGGAGGAAAGTGGATAGTTAGCTATTCATGGGTGGTCCGGTCTATCCAAGAAAGAAGACTTCTGAATGTGCATGAATTTGAAGTCAAAGGAGATGTTGTGACTGGAAGAAATCACCAAGGTCCAAGGCGATCCAGAGAATCCCGGGAAAAGCTCTTCAAGGGCCTACAGGTCTATTGTTGTGAGCCCTTCACCAACATGCCCAAAGATGAGCTGGAGAGGATGCTGCAGCTGTGTGGGGCTTCCGTGGTGAAGGAGCTTCCATCGCTCACCCATGACACAGGTGCTCATCTAGTTGTGATCGTGCAGCCAAGCGCCTGGACAGAAGACAGCAACTGCCCAGATATTGGGCAGCTGTGCAAGGCACGTCTTGTGATGTGGGACTGGGTGTTGGACAGTCTATCCAGCTACCGGTGTCGGGATCTGGATGCCTACCTGGTACAGAATATCACCTGTGACAGTAGTGAGCCACAAGACTCCAATGATTAA'

GAP_PEN_A = -2
GAP_PEN_B = -10

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

SCORING_M = {x: {y: 1 if x==y else -1 if (x in PUR and y in PUR) or (x in PYR and y in PYR) else -2 for y in BASES} for x in BASES}

# Initialisation
for i in range(0,M1):
    D[i,0] = i * GAP_PEN_B
    HA1[i,0], HA2[i,0], HA3[i,0] = float("-inf"), float("-inf"), float("-inf")
for j in range(0,N1):
    D[0,j] = j * GAP_PEN_B
    VA1[0,j], VA2[0,j], VA3[0,j] = float("-inf"), float("-inf"), float("-inf")

print('Recursion...')
# Recursion
for i in range(1,100):
    for j in range(1,100):
        match = SCORING_M[SEQ1[i]][SEQ2[j]]
        Ddiag = D[i-1, j-1] + match
        DtopB = D[i-1, j] + GAP_PEN_B
        DleftB = D[i, j-1] + GAP_PEN_B
        Vtop = VA3[i-1,j-1] + match
        Hleft = HA3[i-1,j-1] + match

        D[i,j] = max(Ddiag, DtopB, DleftB, Vtop, Hleft)

        DtopA = D[i-1, j] + GAP_PEN_A
        VA3topA = VA3[i-1, j] + GAP_PEN_A

        VA1[i,j] = max(DtopA, VA3topA)
        VA2[i,j] = VA1[i-1,j] + GAP_PEN_A
        VA3[i,j] = VA2[i-1,j] + GAP_PEN_A

        DleftA = D[i, j-1] + GAP_PEN_A
        HA3leftA = HA3[i, j-1] + GAP_PEN_A

        HA1[i,j] = max(DleftA, HA3leftA)
        HA2[i,j] = HA1[i, j-1] + GAP_PEN_A
        HA3[i,j] = HA2[i, j-1] + GAP_PEN_A

print(D[90:100,90:100].astype(int))
print(HA3[90:100,90:100].astype(int))
print(VA3[90:100,90:100].astype(int))

# Traceback
print('Tracing back...')
ali1, ali2 = '',''
i, j = 99, 99
state = 'D' if D[i,j] > HA3[i,j] and D[i,j] > VA3[i,j] else 'H3' if HA3[i,j] > VA3[i,j] else 'V3'
score = D[i,j] if state == 'D' else HA3[i,j] if state == 'H3' else VA3[i,j]

while i > 0 and j > 0:
    print('---')
    print(str(i) + '-' + str(j))
    print(state)
    print(score)
    print('-')


    if state == 'D':
        match = SCORING_M[SEQ1[i-1]][SEQ2[j-1]]

        scoreDiag = D[i-1, j-1] + match
        scoreTop = D[i-1, j] + GAP_PEN_B
        scoreLeft = D[i, j-1] + GAP_PEN_B
        scoreH3 = HA3[i, j-1] + match
        scoreV3 = VA3[i-1, j] + match

        print(scoreDiag)
        print(scoreTop)
        print(scoreLeft)
        print(scoreH3)
        print(scoreV3)


        if score == scoreDiag:
            ali1 += SEQ1[i-1]
            ali2 += SEQ2[j-1]
            state = 'D' #redundant
            i -= 1
            j -= 1
            score = D[i,j]
            print('diag')
        elif score == scoreTop:
            ali1 += SEQ1[i-1]
            ali2 += '-'
            state = 'D'
            i -= 1
            score = D[i,j]
            print('top')
        elif score == scoreLeft:
            ali1 += '-'
            ali2 += SEQ2[j-1]
            state = 'D'
            j -= 1
            score = D[i,j]
            print('left')
        elif score == scoreH3:
            ali1 += SEQ1[i-1]
            ali2 += SEQ2[j-1]
            state = 'H3'
            i -= 1
            j -= 1
            score = HA3[i,j]
            print('h3')
        elif score == scoreV3:
            ali1 += SEQ1[i-1]
            ali2 += SEQ2[j-1]
            state = 'V3'
            i -= 1
            j -= 1
            score = VA3[i,j]
            print('v3')
    elif state == 'H1':
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

print(ali1)
print(ali2)

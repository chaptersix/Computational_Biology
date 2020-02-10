import argparse

blosum = {'AA': 4, 'AR': -1, 'AN': -2, 'AD': -2, 'AC': 0, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -1, 'AL': -1, 'AK': -1, 'AM': -1, 'AF': -2, 'AP': -1, 'AS': 1, 'AT':
          0, 'AW': -3, 'AY': -2, 'AV': 0, 'AB': -2, 'AZ': -1, 'AX': 0, 'A*': -4, 'RA': -1, 'RR': 5, 'RN': 0, 'RD': -2, 'RC': -3, 'RQ': 1, 'RE': 0, 'RG': -2, 'RH': 0, 'RI': -3, 'RL': -2, 'RK': 2, 'RM': -1, 'RF': -3, 'RP': -2, 'RS': -1, 'RT': -1, 'RW': -3, 'RY': -2, 'RV': -3, 'RB': -1, 'RZ': 0, 'RX': -1, 'R*': -4, 'NA': -2, 'NR': 0, 'NN': 6, 'ND': 1, 'NC': -3, 'NQ': 0, 'NE': 0, 'NG': 0, 'NH': 1, 'NI': -3, 'NL': -3, 'NK': 0, 'NM': -2, 'NF': -3, 'NP': -2, 'NS': 1, 'NT': 0, 'NW': -4, 'NY': -2, 'NV': -3, 'NB': 3, 'NZ': 0, 'NX': -1, 'N*': -4, 'DA': -2, 'DR': -2, 'DN': 1, 'DD': 6, 'DC': -3, 'DQ': 0, 'DE': 2, 'DG': -1, 'DH': -1, 'DI': -3, 'DL': -4, 'DK': -1, 'DM': -3, 'DF': -3, 'DP': -1, 'DS': 0, 'DT': -1, 'DW': -4, 'DY': -3, 'DV': -3, 'DB': 4, 'DZ': 1, 'DX': -1, 'D*': -4, 'CA': 0, 'CR': -3, 'CN': -3, 'CD': -3, 'CC': 9, 'CQ': -3, 'CE': -4, 'CG': -3, 'CH': -3, 'CI': -1, 'CL': -1, 'CK': -3, 'CM': -1, 'CF': -2, 'CP': -3, 'CS': -1, 'CT': -1, 'CW': -2, 'CY': -2, 'CV': -1, 'CB': -3, 'CZ': -3, 'CX': -2, 'C*': -4, 'QA': -1, 'QR': 1, 'QN': 0, 'QD': 0, 'QC': -3, 'QQ': 5, 'QE': 2, 'QG': -2, 'QH': 0, 'QI': -3, 'QL': -2, 'QK': 1, 'QM': 0, 'QF': -3, 'QP':
          -1, 'QS': 0, 'QT': -1, 'QW': -2, 'QY': -1, 'QV': -2, 'QB': 0, 'QZ': 3, 'QX': -1, 'Q*': -4, 'EA': -1, 'ER': 0, 'EN': 0, 'ED': 2, 'EC': -4, 'EQ': 2, 'EE': 5, 'EG': -2, 'EH': 0, 'EI': -3, 'EL': -3, 'EK': 1, 'EM': -2, 'EF': -3, 'EP': -1, 'ES': 0, 'ET': -1, 'EW': -3, 'EY': -2, 'EV': -2, 'EB': 1, 'EZ': 4, 'EX': -1, 'E*': -4, 'GA': 0, 'GR': -2, 'GN': 0, 'GD': -1, 'GC': -3, 'GQ': -2, 'GE': -2, 'GG': 6, 'GH': -2, 'GI': -4, 'GL': -4, 'GK': -2, 'GM': -3, 'GF': -3, 'GP': -2, 'GS': 0, 'GT': -2, 'GW': -2, 'GY': -3, 'GV': -3, 'GB': -1, 'GZ': -2, 'GX': -1, 'G*': -4, 'HA': -2, 'HR': 0, 'HN': 1, 'HD': -1, 'HC': -3, 'HQ': 0, 'HE': 0, 'HG': -2, 'HH': 8, 'HI': -3,
          'HL': -3, 'HK': -1, 'HM': -2, 'HF': -1, 'HP': -2, 'HS': -1, 'HT': -2, 'HW': -2, 'HY': 2, 'HV': -3, 'HB': 0, 'HZ': 0, 'HX': -1, 'H*': -4, 'IA': -1, 'IR': -3, 'IN':
          -3, 'ID': -3, 'IC': -1, 'IQ': -3, 'IE': -3, 'IG': -4, 'IH': -3, 'II': 4, 'IL': 2, 'IK': -3, 'IM': 1, 'IF': 0, 'IP': -3, 'IS': -2, 'IT': -1, 'IW': -3, 'IY': -1, 'IV': 3, 'IB': -3, 'IZ': -3, 'IX': -1, 'I*': -4, 'LA': -1, 'LR': -2, 'LN': -3, 'LD': -4, 'LC': -1, 'LQ': -2, 'LE': -3, 'LG': -4, 'LH': -3, 'LI': 2, 'LL': 4, 'LK': -2, 'LM': 2, 'LF': 0, 'LP': -3, 'LS': -2, 'LT': -1, 'LW': -2, 'LY': -1, 'LV': 1, 'LB': -4, 'LZ': -3, 'LX': -1, 'L*': -4, 'KA': -1, 'KR': 2, 'KN': 0, 'KD': -1, 'KC': -3, 'KQ': 1, 'KE': 1, 'KG': -2, 'KH': -1, 'KI': -3, 'KL': -2, 'KK': 5, 'KM': -1, 'KF': -3, 'KP': -1, 'KS': 0, 'KT': -1, 'KW': -3, 'KY': -2, 'KV': -2, 'KB': 0, 'KZ': 1, 'KX': -1, 'K*': -4, 'MA': -1, 'MR': -1, 'MN': -2, 'MD': -3, 'MC': -1, 'MQ': 0, 'ME': -2, 'MG': -3, 'MH': -2, 'MI': 1, 'ML': 2, 'MK': -1, 'MM': 5, 'MF': 0, 'MP': -2, 'MS': -1, 'MT': -1, 'MW': -1, 'MY': -1, 'MV': 1, 'MB': -3, 'MZ': -1, 'MX': -1, 'M*': -4, 'FA': -2, 'FR': -3, 'FN': -3, 'FD': -3, 'FC': -2, 'FQ': -3, 'FE': -3, 'FG': -3, 'FH': -1, 'FI': 0, 'FL': 0, 'FK': -3, 'FM': 0, 'FF': 6, 'FP': -4, 'FS': -2, 'FT': -2, 'FW': 1, 'FY': 3, 'FV': -1, 'FB': -3, 'FZ': -3, 'FX': -1, 'F*': -4, 'PA': -1, 'PR': -2, 'PN': -2, 'PD': -1, 'PC': -3, 'PQ': -1, 'PE': -1, 'PG': -2, 'PH': -2, 'PI': -3, 'PL': -3, 'PK': -1, 'PM': -2, 'PF': -4, 'PP': 7, 'PS': -1, 'PT': -1, 'PW': -4, 'PY': -3, 'PV': -2, 'PB': -2, 'PZ': -1, 'PX': -2, 'P*': -4, 'SA': 1, 'SR': -1, 'SN': 1, 'SD': 0, 'SC': -1, 'SQ': 0, 'SE': 0, 'SG': 0, 'SH': -1,
          'SI': -2, 'SL': -2, 'SK': 0, 'SM': -1, 'SF': -2, 'SP': -1, 'SS': 4, 'ST': 1, 'SW': -3, 'SY': -2, 'SV': -2, 'SB': 0, 'SZ': 0, 'SX': 0, 'S*': -4, 'TA': 0, 'TR': -1,
          'TN': 0, 'TD': -1, 'TC': -1, 'TQ': -1, 'TE': -1, 'TG': -2, 'TH': -2, 'TI': -1, 'TL': -1, 'TK': -1, 'TM': -1, 'TF': -2, 'TP': -1, 'TS': 1, 'TT': 5, 'TW': -2, 'TY':
          -2, 'TV': 0, 'TB': -1, 'TZ': -1, 'TX': 0, 'T*': -4, 'WA': -3, 'WR': -3, 'WN': -4, 'WD': -4, 'WC': -2, 'WQ': -2, 'WE': -3, 'WG': -2, 'WH': -2, 'WI': -3, 'WL': -2, 'WK': -3, 'WM': -1, 'WF': 1, 'WP': -4, 'WS': -3, 'WT': -2, 'WW': 11, 'WY': 2, 'WV': -3, 'WB': -4, 'WZ': -3, 'WX': -2, 'W*': -4, 'YA': -2, 'YR': -2, 'YN': -2, 'YD':
          -3, 'YC': -2, 'YQ': -1, 'YE': -2, 'YG': -3, 'YH': 2, 'YI': -1, 'YL': -1, 'YK': -2, 'YM': -1, 'YF': 3, 'YP': -3, 'YS': -2, 'YT': -2, 'YW': 2, 'YY': 7, 'YV': -1, 'YB': -3, 'YZ': -2, 'YX': -1, 'Y*': -4, 'VA': 0, 'VR': -3, 'VN': -3, 'VD': -3, 'VC': -1, 'VQ': -2, 'VE': -2, 'VG': -3, 'VH': -3, 'VI': 3, 'VL': 1, 'VK': -2, 'VM': 1,
          'VF': -1, 'VP': -2, 'VS': -2, 'VT': 0, 'VW': -3, 'VY': -1, 'VV': 4, 'VB': -3, 'VZ': -2, 'VX': -1, 'V*': -4, 'BA': -2, 'BR': -1, 'BN': 3, 'BD': 4, 'BC': -3, 'BQ': 0, 'BE': 1, 'BG': -1, 'BH': 0, 'BI': -3, 'BL': -4, 'BK': 0, 'BM': -3, 'BF': -3, 'BP': -2, 'BS': 0, 'BT': -1, 'BW': -4, 'BY': -3, 'BV': -3, 'BB': 4, 'BZ': 1, 'BX': -1, 'B*': -4, 'ZA': -1, 'ZR': 0, 'ZN': 0, 'ZD': 1, 'ZC': -3, 'ZQ': 3, 'ZE': 4, 'ZG': -2, 'ZH': 0, 'ZI': -3, 'ZL': -3, 'ZK': 1, 'ZM': -1, 'ZF': -3, 'ZP': -1, 'ZS': 0, 'ZT': -1, 'ZW': -3, 'ZY': -2, 'ZV': -2, 'ZB': 1, 'ZZ': 4, 'ZX': -1, 'Z*': -4, 'XA': 0, 'XR': -1, 'XN': -1, 'XD': -1, 'XC': -2, 'XQ': -1, 'XE': -1, 'XG': -1, 'XH': -1, 'XI': -1, 'XL': -1, 'XK': -1, 'XM': -1, 'XF': -1, 'XP': -2, 'XS': 0, 'XT': 0, 'XW': -2, 'XY': -1, 'XV': -1, 'XB': -1, 'XZ': -1, 'XX': -1, 'X*': -4, '*A': -4, '*R': -4, '*N': -4, '*D': -4, '*C': -4, '*Q': -4, '*E': -4, '*G': -4, '*H': -4, '*I': -4, '*L': -4, '*K': -4, '*M': -4, '*F': -4, '*P': -4, '*S': -4, '*T': -4, '*W': -4, '*Y': -4, '*V': -4, '*B': -4, '*Z': -4, '*X': -4, '**': 1}


def open_sequence(filename):
    name = ""
    seq = ""
    with open(filename, 'r') as file:
        name = file.readline()
        seq = file.read().replace("\n", "")
    return(seq)


def prettyPrint(al1, al2, penalties={'m': 1, 'n': -1, "g": -2}):
    score = 0
    bwt = ""
    al1r = ""
    al2r = ""
    for x in range(len(al1)-1, 0, -1):
        if "-" == al1[x] or "-" == al2[x]:
            bwt += " "
            score += penalties['g']
        elif al1[x] == al2[x]:
            bwt += "|"
            score += penalties['m']
        else:
            bwt += "*"
            score += penalties['n']
        al1r += al1[x]
        al2r += al2[x]
    out_str = "Score:{}\n".format(score)
    for x in range(80, len(al1r), 80):
        out_str += "\n{}\n{}\n{}\n".format(al1r[:80], bwt[:80], al2r[:80])
        al1r = al1r[80:]
        bwt = bwt[80:]
        al2r = al2r[80:]
    out_str += "\n{}\n{}\n{}\n".format(al1r, bwt, al2r)
    return(out_str)


def zeroMatrix(n, m):
    matrix = []
    for x in range(n):
        matrix.append([0]*m)
    return(matrix)


def NW(seq1, seq2, penalties={'m': 5, 'n': -1, "g": -2}):
    n = len(seq1) + 1
    m = len(seq2) + 1
    matrix = zeroMatrix(n, m)

    for r in range(len(matrix)):
        matrix[r][0] = r * penalties["g"]
    for c in range(len(matrix[0])):
        matrix[0][c] = c * penalties["g"]

    def score_match(a, b):
        return(blosum[a+b])

    for r in range(1, n):
        for c in range(1, m):

            Diagonal = matrix[r-1][c-1] + score_match(seq1[r-1], seq2[c-1])
            Insert = matrix[r][c-1] + penalties['g']
            Del = matrix[r-1][c] + penalties['g']

            matrix[r][c] = max(Diagonal, Insert, Del)

    aligned_seq1 = ""
    aligned_seq2 = ""
    i = n - 1
    j = m - 1
    while i > 0 and j > 0:
        score_current = matrix[i][j]
        score_up = matrix[i-1][j]
        score_left = matrix[i][j-1]
        score_diag = matrix[i-1][j-1]

        if score_current == score_diag + score_match(seq1[i-1], seq2[j-1]):
            seq1_char = seq1[i-1]
            seq2_char = seq2[j-1]
            i = i-1
            j = j-1
        elif score_current == score_up + penalties['g']:
            seq1_char = seq1[i-1]
            seq2_char = "-"
            i += -1
        elif score_current == score_left + penalties['g']:
            seq1_char = '-'
            seq2_char = seq2[j-1]
            j += -1
        aligned_seq1 += seq1_char
        aligned_seq2 += seq2_char

    return(prettyPrint(aligned_seq1, aligned_seq2, penalties))


def SW(seq1, seq2, penalties={'m': 8, 'n': -1, "g": -10}):
    n = len(seq1) + 1
    m = len(seq2) + 1
    matrix = zeroMatrix(n, m)

    def maxRow(i, j):
        maxInt = matrix[i][j]
        while j > 0:
            if matrix[i][j] > maxInt:
                maxInt = matrix[i][j]
            j -= 1
        return(maxInt)

    def maxCol(i, j):
        maxInt = matrix[i][j]
        while j > 0:
            if matrix[i][j] > maxInt:
                maxInt = matrix[i][j]
            j -= 1
        return(maxInt)
    # Fill
    max_location = (0, [0, 0])
    for r in range(1, n):
        for c in range(1, m):
            match = matrix[r-1][c-1] + blosum[seq1[r-1] +
                                              seq2[c-1]] if seq1[r-1] == seq2[c-1] else 0
            insert = maxRow(r, c) - penalties['g']
            delete = maxCol(r, c) - penalties['g']
            maxed = max(match, insert, delete, 0)
            if maxed > max_location[0]:
                max_location = (maxed, [r, c])
            matrix[r][c] = maxed

    i = max_location[1][0]
    j = max_location[1][1]
    aligned_seq1 = ""
    aligned_seq2 = ""
    while i > 0 or j > 0:
        dia = ("d", matrix[i-1][j-1])
        top = ("t", matrix[i-1][j])
        left = ("l", matrix[i][j-1])
        maxed = max(dia, top, left, key=lambda x: x[1])

        if maxed[0] == 'd':
            seq1_char = seq1[i-1]
            seq2_char = seq2[j-1]
            i -= 1
            j -= 1
        elif maxed[0] == 't':
            seq1_char = seq1[i-1]
            seq2_char = '-'
            i -= 1
        elif maxed[0] == 'l':
            seq1_char = '-'
            seq2_char = seq2[j-1]
            j -= 1
        aligned_seq1 += seq1_char
        aligned_seq2 += seq2_char
        if matrix[i-1][j-1] == 0:
            break
    return(prettyPrint(aligned_seq1, aligned_seq2, penalties))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("sequence_one_filename",
                        help="relative path to the file containing sequence one")
    parser.add_argument("sequence_two_filename",
                        help="relative path to the file containing sequence two")
    parser.add_argument("output_filename",
                        help="file that you would like the output written to")
    args = parser.parse_args()

    seq1 = open_sequence(args.sequence_one_filename)
    seq2 = open_sequence(args.sequence_two_filename)
    nw_solution = NW(seq1, seq2)
    sw_solution = SW(seq1, seq2)
    with open(args.output_filename, "w") as file:
        file.write("Global Alignment\n{}\nLocal Alignment\n{}".format(
            nw_solution, sw_solution))

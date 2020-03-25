import csv
import math
import argparse


def read_ss(protein_name):
    ss = []
    with open('ss\\\\'+protein_name+'.ss', 'r') as file:
        ss = list(file.readlines()[1].replace('\n', ''))
    return(ss)


def read_total_probs():
    p = {}
    with open('total_probs.txt', 'r') as file:
        p = eval(file.read())
    return(p)


def read_means_variances():
    p = {}
    with open('means_variances.txt', 'r') as file:
        p = eval(file.read())
    return(p)


def read_pssm_file(protein_name):
    ss = read_ss(protein_name)
    matrix = []
    with open('pssm\\\\'+protein_name+'.pssm', 'r') as file:
        reader = csv.reader(file, delimiter=' ')
        for num, row in enumerate(reader):
            if num > 2:
                m_row = list(filter(lambda x: x != "", row))[2:22]
                if len(m_row) <= 0:
                    break
                m_tuple = (m_row, ss[num-3])
                matrix.append(m_tuple)
    return(matrix)


def sliding_window(matrix):
    windowed_matrix = []

    def make_list(matrix_rows):
        re = []
        for l in matrix_rows:
            re = re + l
        return(re)

    for num, row in enumerate(matrix):
        if num > 1 and num < len(matrix) - 2:
            nrow = (make_list([x[0] for x in matrix[num-2:num+3]]),
                    row[1])  # average case
        if num < 2:
            l = make_list([x[0] for x in matrix[num:num+3]])
            negs = [-1 for _ in range(100-len(l))]
            nrow = (negs+l, row[1])
        if num >= len(matrix)-2:
            l = make_list([x[0] for x in matrix[num-2:]])
            negs = [-1 for _ in range(100-len(l))]
            nrow = (l+negs, row[1])

        windowed_matrix.append(nrow)
    return(windowed_matrix)


def process_protein(protein_name):
    m = read_pssm_file(protein_name)
    wm = sliding_window(m)
    return(wm)


def gnb(xi, class_mean, class_variance, class_prob):
    bottom = math.sqrt(2*math.pi*class_variance)
    r = math.pow((xi-class_mean)/math.sqrt(class_variance), 2)
    r = -0.5 * r
    top = math.exp(r)
    return(top/bottom)


def classify(protein_name):
    m = process_protein(protein_name)
    t_probs = read_total_probs()
    mv = read_means_variances()
    cls = ["H", "E", "C"]
    output = ""
    for row in m:
        col_probs = {}
        for i in range(100):
            for cl in cls:
                if cl not in col_probs:
                    col_probs[cl] = gnb(
                        int(row[0][i]), mv[i][cl][0], mv[i][cl][1], t_probs[cl])
                else:
                    col_probs[cl] *= gnb(int(row[0][i]), mv[i]
                                         [cl][0], mv[i][cl][1], t_probs[cl])
        for k in cls:
            col_probs[k] *= t_probs[k]
        high_n = list(col_probs.items())[0][0]
        high_v = list(col_probs.items())[0][1]
        for k, v in col_probs.items():
            if v > high_v:
                high_n = k
                high_v = v
        output += high_n
    return(output)


def test():
    test_proteins = ['1vp6', '1fvg', '1dmg', '1jkx', '1bdo', '1ag6', '1gz2', '3bor', '1m8a', '1brf', '1ej0', '1aba', '1o1z', '1kw4', '1ktg', '3dqg', '1ihz', '1cke',
                     '1jl1', '1dix', '1wjx', '2phy', '1jo8', '1kqr', '1mk0', '1vhu', '1fx2', '1htw', '1beh', '1aap', '1kq6', '1h4x', '1c9o', '1c44', '1tzv', '1qf9', '1bsg', '1smx']

    total_correct = 0
    count = 0
    for protein_name in test_proteins:
        expected = list(read_ss(protein_name))
        actual = list(classify(protein_name))
        for x in range(len(actual)):
            if expected[x] == actual[x]:
                total_correct += 1
            count += 1
    return(total_correct/count)


def read_fasta_protein_name(file_path):
    name = ""
    with open(file_path, "r") as file:
        name = file.readline()[1:-1]
    return(name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the project classifier')
    parser.add_argument('fasta_file', metavar='Fasta_file', type=str,
                        help='path to the fasta file')
    args = parser.parse_args()
    name = read_fasta_protein_name(args.fasta_file)
    print(classify(name))

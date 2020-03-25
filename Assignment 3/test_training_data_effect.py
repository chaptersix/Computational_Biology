import math
import csv
import glob
import random
from tqdm import tqdm


def split_data(r):
    files = glob.glob('fasta/*')
    files = [x.replace("fasta\\", "").replace(".fasta", "") for x in files]
    r.shuffle(files)

    training_proteins = files[:112]
    testing_proteins = files[112:]
    return((training_proteins, testing_proteins))


def read_ss(protein_name):
    ss = []
    with open('ss\\\\'+protein_name+'.ss', 'r') as file:
        ss = list(file.readlines()[1].replace('\n', ''))
    return(ss)


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


def generate_feature_matrix(training_proteins):
    fm = []
    for name in training_proteins:
        wm = process_protein(name)
        for r in wm:
            fm.append(r)
    return(fm)


def col_mean_variance(matrix, col_number, cl):
    mean = 0
    mean_count = 0
    for row in matrix:
        el = row[1]
        if el == cl:
            mean += int(row[0][col_number])
            mean_count += 1
    mean = mean / mean_count
    variance = 0
    variance_count = 0
    for row in matrix:
        el = row[1]
        if el == cl:
            variance += math.pow(int(row[0][col_number])-mean, 2)
            variance_count += 1
    variance = variance / variance_count
    return(mean, variance)


def total_probs(matrix, cls):
    cls_probs = {}
    l = [x[1] for x in matrix]
    for cl in cls:
        nums = l.count(cl)
        avg = nums/len(l)
        cls_probs[cl] = avg
    return(cls_probs)


def generate_means_variances(matrix):
    cls = ["H", "E", "C"]
    # cls_prob = total_probs(matrix, cls)
    col_values = []
    for i in range(100):
        cls_col_prob = {}
        for cl in cls:
            cls_col_prob[cl] = col_mean_variance(matrix, i, cl)

        col_values.append(cls_col_prob)

    return(col_values)


def gnb(xi, class_mean, class_variance, class_prob):
    bottom = math.sqrt(2*math.pi*class_variance)
    r = math.pow((xi-class_mean)/math.sqrt(class_variance), 2)
    r = -0.5 * r
    top = math.exp(r)
    return(top/bottom)


def classify(protein_name, mv, t_probs):
    m = process_protein(protein_name)
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


def train_test(random_class):
    train, test = split_data(random_class)
    m = generate_feature_matrix(train)
    means_variances = generate_means_variances(m)
    total_prob = total_probs(m, ["H", "E", "C"])

    total_correct = 0
    count = 0
    for protein_name in test:
        expected = list(read_ss(protein_name))
        actual = list(classify(protein_name, means_variances, total_prob))
        for x in range(len(actual)):
            if expected[x] == actual[x]:
                total_correct += 1
            count += 1
    return(total_correct/count)


iters = []
r = random.Random(4)
for x in tqdm(range(100)):
    iters.append(train_test(r))
print(iters)

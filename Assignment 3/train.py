import csv
import math
train_names = ['2cua', '1nps', '2arc', '1jvw', '1im5', '1jbe', '1ctf', '1a3a', '1h0p', '1h2e', '1gbs', '1xkr', '1i58', '1whi', '1avs', '1svy', '1jwq', '1rw7', '1d0q', '1nb9', '1ny1', '1dlw', '1xdz', '1w0h', '2tps', '1jo0', '1bkr', '1chd', '1ek0', '1vjk', '1cjw', '1g2r', '1tqh', '1pko', '1mug', '1dqg', '1tqg', '1pch', '1hxn', '1g9o', '1lpy', '1jfx', '1rw1', '1i1j', '1c52', '1fcy', '1d1q', '1m4j', '1nrv', '1xff', '1dbx', '1wkc', '1j3a', '1fk5', '1tif',
               '1vmb', '1i4j', '1fl0', '1i1n', '1beb', '1jbk', '1jfu', '1h98', '1iwd', '1lm4', '1vfy', '1p90', '1a6m', '1gzc', '1hh8', '1iib', '1ej8', '2hs1', '1czn', '1fvk', '1k7c', '1qjp', '1ku3', '2vxn', '1aoe', '1jyh', '1fqt', '1jos', '5ptp', '1i5g', '1atl', '1r26', '1ne2', '1k7j', '1a70', '1hfc', '1f6b', '1cc8', '1atz', '1d4o', '1hdo', '1gmi', '1lo7', '1k6k', '1gmx', '1kid', '1fna', '1i71', '1dsx', '1guu', '1ql0', '1roa', '1ryb', '1cxy', '1t8k', '2mhr', '1eaz']


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


m = generate_feature_matrix(train_names)
means_variances = generate_means_variances(m)
total_probs = total_probs(m, ["H", "E", "C"])

# with open("means_variances.txt", "w") as file:
#     file.write(str(means_variances))


# with open("total_probs.txt", "w") as file:
#     file.write(str(total_probs))

with open('feat.txt', "w") as file:
    file.write(str(m))

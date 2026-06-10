import os

import pandas as pd


PYTHONDIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BARCODE_XLSX = os.path.join(PYTHONDIR, "Barcode.xlsx")
W1 = "GAGTGATTGCTTGTGACGCCTT"
BC2_NUMBER = 384
BARCODE_DIGITS = 5


_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "n": "n",
}


def reverse_complement(seq):
    return "".join(_complement.get(base, "N") for base in reversed(seq.upper()))


def hamming_distance(seq1, seq2):
    return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))


def find_w1_loc(read, w1=W1, threshold=5):
    length_w1 = len(w1)
    start_loc = 8
    tests = 4
    closest_loc = start_loc
    diff = length_w1

    for i in range(tests):
        candidate = read[start_loc + i:start_loc + i + length_w1]
        if len(candidate) != length_w1:
            continue

        diff_temp = hamming_distance(candidate, w1)
        if diff > diff_temp:
            closest_loc = start_loc + i
            diff = diff_temp

        if diff < threshold:
            return closest_loc

    return 0


def closest_index(seq, ref_list, mismatch_allow=1):
    closest_index_list = []
    closest_ham_dist = len(seq)

    for ref in ref_list:
        ham_dist_temp = hamming_distance(seq, ref)
        if ham_dist_temp < closest_ham_dist:
            closest_ham_dist = ham_dist_temp
            closest_index_list = [ref]
        elif ham_dist_temp == closest_ham_dist:
            closest_index_list.append(ref)

    if closest_ham_dist <= mismatch_allow:
        return closest_index_list

    return []


def load_barcode_reference(barcode_xlsx=DEFAULT_BARCODE_XLSX, bc1_number=96, bc2_number=BC2_NUMBER):
    barcode_temp = pd.read_excel(barcode_xlsx)

    bc1_temp = barcode_temp["bc1"].values
    bc2_temp = barcode_temp["bc2"].values

    bc2 = {}
    for i in range(bc2_number):
        bc2[str(bc2_temp[i][31:39]).upper()] = i

    bc1 = {}
    for i in range(bc1_number):
        loc = bc1_temp[i].find("AGATCGGAAGAGCGTCGTGTAGGGAAAGAG")
        bc1[str(bc1_temp[i][23:loc - 1]).upper()] = i

    return bc1, bc2


bc1, bc2 = load_barcode_reference()
bc1_list = list(bc1.keys())
bc2_list = list(bc2.keys())


def find_barcode(read):
    read = read.strip().upper()

    if len(read) <= 80:
        return "Len"

    loc = find_w1_loc(read)
    if loc == 0:
        return "W1"

    bc1_temp = reverse_complement(read[:loc])
    bc2_temp = reverse_complement(read[loc + 22:loc + 30])

    if bc1_temp in bc1_list:
        read_bc1 = bc1[bc1_temp]
    else:
        bc1_cl = closest_index(bc1_temp, bc1_list)
        if len(bc1_cl) == 1:
            read_bc1 = bc1[bc1_cl[0]]
        else:
            return "BC"

    if bc2_temp in bc2_list:
        read_bc2 = bc2[bc2_temp]
    else:
        bc2_cl = closest_index(bc2_temp, bc2_list)
        if len(bc2_cl) == 1:
            read_bc2 = bc2[bc2_cl[0]]
        else:
            return "BC"

    return read_bc1 * BC2_NUMBER + read_bc2


def barcode_tag_func(read, sample_id):
    barcode = find_barcode(read)

    if isinstance(barcode, int):
        return f"Cell{sample_id}{barcode:0{BARCODE_DIGITS}d}", "Valid"

    return barcode, f"Error:{barcode}"

def make_tag_func(sample_id):
    def tag_func(read):
        return barcode_tag_func(read, sample_id)
    return tag_func
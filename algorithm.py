# import glob


propensity_file = 'data/structure_propensities.txt'
sequence_file = 'data/test_sequence.fasta'


def read_sequence(file):
    file = open(file, 'r')
    aa_sequence = ''
    for line in file:
        if line.startswith('>'):
            pass
        else:
            line = line.replace('\n', '').replace('\r', '')
            aa_sequence += line
    return aa_sequence


# creates dictionary and adds AS (key) + desired colum (value) into it
def read_propensities_into_dic(file, column):
    file = open(file, 'r')
    dic_temp = {}
    word_matrix = []
    word = ''
    for line in file:
        if line.startswith('-'):
            pass
        else:
            word_list = []
            for char in line:
                if char != ' ' and char != '\n' and char != '\r':  # \r for linux
                    word = word + char
                elif (char == ' ' or char == '\n' or char == '\r') and word != '':
                    # last value missing if it's last line of file. fixed by adding another '--' line.
                    word_list.append(word)
                    word = ''
                else:
                    pass
            word_matrix.append(word_list)
    for element in word_matrix:
        dic_temp[element[0]] = float(element[column])
    return dic_temp


print('Sequence: ', read_sequence(sequence_file))
sequence = read_sequence(sequence_file)

helix_dic = read_propensities_into_dic(propensity_file, 1)  # opens file 3 times -> should be changed
sheet_dic = read_propensities_into_dic(propensity_file, 2)
turn_dic = read_propensities_into_dic(propensity_file, 3)
print('Helix: ', helix_dic)
print('Sheet: ', sheet_dic)
print('Turn: ', turn_dic)


# creates list of secondary structure values for sequence
def map_dic_to_seq(sequence, dictionary):
    value_list = []
    for aa in sequence:
        value_list.append(dictionary[aa])
    return value_list


helix_value_list = (map_dic_to_seq(sequence, helix_dic))
sheet_value_list = (map_dic_to_seq(sequence, sheet_dic))
turn_value_list = (map_dic_to_seq(sequence, turn_dic))
print('Helix_value_list: ', helix_value_list)
print('Sheet_value_list: ', sheet_value_list)
print('Turn_value_list: ', turn_value_list)


def find_extendindex(motif, strseq):
    index_list = []
    stpoint = 0
    position = 0
    while position >= 0:
        position = strseq.find(motif, stpoint, len(strseq))
        if position >= 0:
            index_list.append(position)
        stpoint = position + len(motif)
    return index_list


def find_nucleations(value_list, secondary_structure):
    if secondary_structure == 'H':
        window_size = 6
        amount_for_nuc = 4
        nucleation_threshold = 1.03
    elif secondary_structure == 'E':
        window_size = 5
        amount_for_nuc = 3
        nucleation_threshold = 1.00
    # elif secondary_structure == 'T':
    #    pass
    sec_struct_list = ['-'] * len(value_list)
    for aa in range(len(value_list) - window_size):
        window = []
        for win_index in range(aa, aa + window_size):  # index error somewhere in here??? line 75 & 79
            window.append(value_list[win_index])
        result = list(filter(lambda x: x >= nucleation_threshold, window))  # look at documentation if only > or >=
        if len(result) >= amount_for_nuc:
            for i in range(aa, aa + window_size):
                sec_struct_list[i] = secondary_structure   # may not be quite correct according to original paper
    stringseq = ''.join(sec_struct_list)
    # print('HIER!:', stringseq)
    fwd =secondary_structure * 3 + '-'
    rev ='-' + secondary_structure * 3
    # print(fwd, rev)
    print('Liste von:', secondary_structure)
    print('fwd', find_extendindex(fwd, stringseq))
    print('rev', find_extendindex(rev, stringseq))
    checkseq = stringseq
    loopcheck = True
    while stringseq != checkseq or loopcheck is True:
        fwdindex = find_extendindex(fwd, stringseq)
        revindex = find_extendindex(rev, stringseq)
        # print('hier!', fwdindex)
        #print('hier!', revindex)
        checkseq = stringseq
        loopcheck = False
        # erweiterung
        #print('value list', value_list)
        for i in fwdindex:  # fwd
            fwd_threshold = 0
            for j in range(i, i + 4):
                fwd_threshold = fwd_threshold + value_list[j]
            fwd_threshold = fwd_threshold / 4
            #print('fwd_thresh', fwd_threshold)
            if fwd_threshold > 1.00:
                stringseq = stringseq[:i + 1] + secondary_structure * 4 + stringseq[i + 5:]

        for i in revindex:  # rev
            rev_threshold = 0
            for j in range(i, i + 4):
                rev_threshold = rev_threshold + value_list[j]
            rev_threshold = rev_threshold / 4
            if rev_threshold > 1.00:
                stringseq = stringseq[:i + 1] + secondary_structure * 4 + stringseq[i + 5:]
        #print(checkseq)
        #print(sec_struct_list)
        sec_struct_list = list(stringseq)   # might not be needed


    return sec_struct_list

def extend_nucleations():  # might not be needed/ included in find nucleations/ called by find nucleations??
    pass




helix_nuc_list = find_nucleations(helix_value_list, 'H')
sheet_nuc_list = find_nucleations(sheet_value_list, 'E')
print('Helix_nuc_list: ', helix_nuc_list)
print('Sheet_nuc_list: ', sheet_nuc_list)


# http://www.biogem.org/tool/chou-fasman/index.php
# comparison to webtool with implemented cf
test = 'HHHHHHHHHHHHHHHHHHH                   HHHHHHHHHHHHHHHHHHHHHH       HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH          HHHHHHHHHHHHHHHHHHHH        HHHHHHH          '
test2 = []
for i in test:
    test2.append(i)
print('orig')
print(test2)
print('pred')
print(helix_nuc_list)

testsh = 'EEEEEEE     EEEEEEEEEEEEEEEEEEEE                               EEEEEE       EEE                          EEEEEEEEEEE      EEEEEEEEEEEE              EEEE       EEEEEEEE      '
test2sh = []
for i in testsh:
    test2sh.append(i)
print('orig')
print(test2sh)
print('pred')
print(sheet_nuc_list)



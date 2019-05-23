# import glob


propensity_file = 'data/structure_propensities.txt'
sequence_file = 'data/test_sequence.fasta'

def read_sequence(file):
    file = open(file, 'r')
    sequence = ''
    for line in file:
        if line.startswith('>'):
            pass
        else:
            line = line.replace('\n','').replace('\r','')
            sequence += line
    return sequence

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

helix_dic = read_propensities_into_dic(propensity_file, 1) # opens file 3 times -> should be changed
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

def find_nucleations(value_list,secondary_structure):
    if secondary_structure == 'H':
        window_size = 6
        amount_for_nuc = 4
        nucleation_threshold = 1.03
    elif secondary_structure == 'E':
        window_size = 5
        amount_for_nuc = 3
        nucleation_threshold = 1.00
    #elif secondary_structure == 'T':
    #    pass
    sec_struct_list = ['-'] * len(value_list)
    for aa in range(0, len(value_list) - window_size):
        window = []
        for i in range(aa, aa + window_size):
            window.append(value_list[i])
        result = list(filter(lambda x: x >= nucleation_threshold, window)) # look at documentation if only > or >=
        if len(result) >= amount_for_nuc:
            for i in range(aa, aa + window_size):
                sec_struct_list[i] = secondary_structure
    return sec_struct_list

helix_nuc_list = find_nucleations(helix_value_list, 'H')
sheet_nuc_list = find_nucleations(sheet_value_list, 'E')
print('Helix_nuc_list: ', helix_nuc_list)
print('Sheet_nuc_list: ', sheet_nuc_list)

def extend_nucleations():
    pass


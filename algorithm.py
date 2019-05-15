
# reading Structure_propensities.txt

filename = 'Structure_propensities.txt'


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


helix_dic = read_propensities_into_dic(filename, 1)
sheet_dic = read_propensities_into_dic(filename, 2)
turn_dic = read_propensities_into_dic(filename, 3)
print('Helix', helix_dic)
print('Sheet', sheet_dic)
print('Turn', turn_dic)


# creates list of secondary structure values for sequence
def map_dic_to_seq(sequence, dictionary):
    value_list = []
    for key in sequence:
        value_list.append(dictionary[sequence] = dictionary.values())  # doesn't work yet
    return value_list

print(map_dic_to_seq('ARD', helix_dic))
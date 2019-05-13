
# reading Structure_propensities.txt

filename = 'Structure_propensities.txt'

def read_propensities_into_dic(filename, column): # creates dictionary and adds AS (key) + desired colum (value) into it
    file = open(filename, 'r')
    dic_temp = {}
    word_matrix = []
    word = ''
    for line in file:
        if line.startswith('-'):
            pass
        else:
            word_list = []
            for char in line:
                if char != ' ' and char != '\n' and  char != '\r': # \r for linux
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

print(read_propensities_into_dic(filename,1))

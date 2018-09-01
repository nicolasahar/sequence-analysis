## Section 1.3 - Finding mutations using a reference genome ##

import os
os.chdir("/Users/nicolasahar/Desktop/ECs/Archive/2015-2016/Computing for Medicine (Apr 2016)/Phase 3/Seminar 3/Data")

from collections import Counter

def find_mutations(file_name):
    answer = [file_name, 0, []]

    r = open("EM_079517.fasta", "rU")
    r.readline()
    reference = r.readline()

    f = open(file_name, "rU")
    f.readline()
    current = "".join(line.strip() for line in f)

    for char in range(len(current)):
        if current[char] == "N" or current[char] == "\n":
            pass

        elif current[char] != reference[char]:
            answer[1] += 1
            answer[2].append([char+1, reference[char], current[char]])

    return len(reference), answer

print(find_mutations("EBOV_20236_MinION_GUI_Boke_2015-06-20.genome.fasta"))

def final_mutations():
    file_list = os.listdir()

    answer_list = [] # each entry is a list with the following format: <file name>, <total number of mutations>, <list of all mutations>
    for file in file_list:
        if ".genome.fasta" in file:
            answer_list.append(find_mutations(file))

    #finding the most common mutation
    genome_list = [x[2] for x in answer_list]

    mutation_list = []
    for genome in genome_list:
        for mutation in genome:
            mutation_list.append(mutation[0])

    data = Counter(mutation_list)
    print(data.most_common())

    return answer_list

## Section 1.4 - Clustering samples ##

def find_common_mutations(genome1, genome2):
    answer = [[genome1, genome2], 0, []]

    g1 = open(genome1, "rU")
    g1.readline()
    gen1 = "".join(line.strip() for line in g1)

    g2 = open(genome2, "rU")
    g2.readline()
    gen2 = "".join(line.strip() for line in g2)

    for char in range(len(gen1)):
        if gen1[char] == "N" or gen2[char] == "N":
            pass

        elif gen1[char] != gen2[char]:
            answer[1] += 1
            answer[2].append([char+1, gen1[char], gen2[char]])

    return answer

# General function that finds all mutations between all sequences
def final_common_mutations():
    file_list = os.listdir()

    genome_list = []

    answer_list = [] # each entry is a list with the following format: <file name1>, <file name2>, <total number of mutations>, <list of all mutations>
    for file in file_list:
        if ".genome.fasta" in file:
            genome_list.append(file)

    for g1 in range(len(genome_list)):
        for g2 in range(g1+1,len(genome_list)):
            answer_list.append(find_common_mutations(genome_list[g1], genome_list[g2]))

    return answer_list

# Function to answer: What genome is closest to sample EBOV_EM_COY_2015_017498_MinION_GUI_Forecariah_2015-06-08?
def question_1():
    file_list = os.listdir()

    genome_list = []

    answer_list = [] # each entry is a list with the following format: <file name1>, <file name2>, <total number of mutations>, <list of all mutations>
    for file in file_list:
        if ".genome.fasta" in file:
            genome_list.append(file)

    g1 = "EBOV_EM_COY_2015_017498_MinION_GUI_Forecariah_2015-06-08.genome.fasta"

    for g2 in genome_list:

        ans = find_common_mutations(g1, g2)
        answer_list.append(ans)

        print(ans[1],ans)

    return answer_list

#print(question_1())

# Function to answer: How often do a pair of samples have the exact same sequence (no positions are different). What might this tell you about these cases?
def same_sequence():
    file_list = os.listdir()

    genome_list = []

    same = 0
    total = 0

    for file in file_list:
        if ".genome.fasta" in file:
            genome_list.append(file)

    for g1 in range(len(genome_list)):
        for g2 in range(g1+1,len(genome_list)):

            ans = find_common_mutations(genome_list[g1], genome_list[g2])

            #print(ans[1], ans)

            total += 1

            if ans[1] == 0:
                same +=1

    return same, total, same/total

#print(same_sequence())
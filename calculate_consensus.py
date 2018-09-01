## Section 1.2 Reconstruct a genome from sequencing reads ##

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import glob

import os
os.chdir("/Users/nicolasahar/Desktop/ECs/Archive/2015-2016/Computing for Medicine (Apr 2016)/Phase 3/Seminar 3/Data")

def calculate_consensus_from_reads(input_filename):
    
    # Open the file containing the input reads
    handle = open(input_filename, "rU")

    # Iterate over the input reads and save their sequence in a list
    read_sequences = []
    for record in SeqIO.parse(handle, "fasta"):
        read_sequences.append(str(record.seq))

    # Initialize the genome sequence as a list of Ns (unknown nucleotide)
    genome_length = len(read_sequences[0])
    genome_sequence = ['N'] * genome_length

    # Write code in the next section to calculate the nucleotide at each position
    # of the genome by finding the most frequent base at that position in the input reads
    # genome_sequence[i] = get_most_frequent_base(read_sequences, i)

    for i in range(len(read_sequences[0])):
        count_list = [0, 0, 0, 0, 0] #keeps count of "N", "A", "T", "G", and "C" respectively
        char_list = ["N", "A", "T", "G", "C"]

        for sequence in read_sequences:
            n = sequence[i]

            if n == "N":
                count_list[0] += 1

            elif n == "A":
                count_list[1] += 1

            elif n == "T":
                count_list[2] += 1

            elif n == "G":
                count_list[3] += 1

            elif n == "C":
                count_list[4] += 1

        max_index = count_list.index(max(count_list))
        genome_sequence[i] = char_list[max_index]

    #return genome_sequence

    # Write the result to disk

    # Use the name of the reads file as the root of the output file name
    out_name = input_filename.replace(".reads.fasta", ".genome.fasta")

    # Open the output file
    output_handle = open(out_name, "w")

    # Make a SeqRecord object for the consensus sequence and write to the output file.
    out_record = SeqRecord(Seq(''.join(genome_sequence)),
                           id=out_name,
                           description="")

    SeqIO.write(out_record, output_handle, "fasta")

# Find all files with the name "*.reads.fasta" in the current working directory and run calculate_consensus_from_reads on those files

def final_consensus():
    file_list = os.listdir()

    for file in file_list:
        if ".reads.fasta" in file:
            calculate_consensus_from_reads(file)
print(final_consensus())



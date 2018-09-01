## Section 1.5 - Building a phylogenetic tree ##

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import glob

os.chdir("/Users/nicolasahar/Desktop/ECs/Archive/2015-2016/Computing for Medicine (Apr 2016)/Phase 3/Seminar 3/Data")

output_handle = open("all_genomes.fasta", "w")

for filename in glob.glob("*.genome.fasta"):
    input_handle = open(filename, "rU")
    record = SeqIO.parse(input_handle, "fasta")
    SeqIO.write(record, output_handle, "fasta")

output_handle.close()

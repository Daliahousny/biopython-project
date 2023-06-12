from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os.path
import sys
import getopt


def gc_percentage(seq):
    seq = Seq(seq)
    gc_percent = gc_fraction(seq)
    print("gc = ", gc_percent * 100)


def dna_transcription(dna):
    dna = Seq(dna)
    transcription = dna.transcribe()
    print("Transcribed DNA: ", transcription)


def rev_complement(dna):
    dna = Seq(dna)
    rev_comp = dna.reverse_complement()
    print("Reversed Complement: ", rev_comp)


def calc_nbases(seq):
    there_is_n = False
    nbase_counter = 0
    for i in seq:
        if i == 'n' or i == "N":
            there_is_n = True
            nbase_counter += 1
    if there_is_n == False:
        print("this sequence doesnt contain any n bases")
    print("Count Nbases = ", nbase_counter)


def is_valid(seq, type):
    dna = 'ATCG'
    rna = 'AUCG'
    amino = 'ARNDCEQGHILKMFPSTWYVUO'
    flag = True

    if type.upper() not in ['DNA', 'RNA', 'PROTEIN']:
        print("Incorrect Type")

    else:
        for i in seq.upper():
            if type.upper() == 'DNA':
                if i not in dna:
                    flag = False
                    break

            elif type.upper() == 'RNA':
                if i not in rna:
                    flag = False
                    break

            elif type.upper() == 'PROTEIN':
                if i not in amino:
                    flag = False
                    break
    print(flag)


def filter_nbases(seq):
    seq_with_noNBases = ''
    for i in seq:
        if i not in ['N', 'n']:
            seq_with_noNBases += i
    print("Sequence with no N bases: ", seq_with_noNBases)


def output_alignment(alignments, output):
    f = open(output, 'w')
    for alignment in alignments:
        not_Format_Alignment = str(alignment)
        f.write(not_Format_Alignment)
        f.write('\n')
        formate_Alignment = str(format_alignment(*alignment))
        f.write(formate_Alignment)
        f.write('\n')
    f.close()


def seq_alignment(seq1, seq2, output=" "):
    seq1 = Seq(seq1)
    seq2 = Seq(seq2)

    alignments = pairwise2.align.globalxx(seq1, seq2)  # global alignment
    if output == " ":
        for alignment in alignments:
            print(alignment)
            print(format_alignment(*alignment))
    else:
        output_alignment(alignments, output)
        print("the output of this alignment in the txt file----> ", output)


def seq_alignment_files(file1, file2, output=""):
    try:
        seq1 = SeqIO.read(file1, "fasta")
        seq2 = SeqIO.read(file2, "fasta")
    except:
        print('Please Enter a valid File name')
        return
    alignments = pairwise2.align.globalxx(seq1, seq2)  # global alignment
    if output == "":
        for alignment in alignments:
            print(alignment)
            print(format_alignment(*alignment))
    else:
        output_alignment(alignments, output)
        print("the output of this alignment in the txt file----> ", output)


def write_blast(blast_record, output):
    f = open(output, 'w')
    for description in blast_record.descriptions:
        f.write("***************************Description***********************************\n")
        f.write("Title " + str(description.title) + '\n')
        f.write("Score " + str(description.score) + '\n')
        f.write("e " + str(description.e) + '\n')
        f.write("Number of Alignments " + str(description.num_alignments) + '\n')
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            f.write("*************************Alignments**********************\n")
            f.write("Title " + str(alignment.title) + '\n')
            f.write("Length " + str(alignment.length) + '\n')
            f.write("****************************hsp*****************************\n")
            f.write("score " + str(hsp.score) + '\n')
            f.write("bits " + str(hsp.bits) + '\n')
            f.write("expect " + str(hsp.expect) + '\n')
            f.write("alignments " + str(hsp.num_alignments) + '\n')
            f.write("identities " + str(hsp.identities) + '\n')
            f.write("positives " + str(hsp.positives) + '\n')
            f.write("gaps " + str(hsp.gaps) + '\n')
            f.write("strand " + str(hsp.strand) + '\n')
            f.write("frame " + str(hsp.frame) + '\n')
            f.write("query " + str(hsp.query) + '\n')
            f.write("query_start " + str(hsp.query_start) + '\n')
            f.write("match " + str(hsp.match) + '\n')
            f.write("sbjct " + str(hsp.sbjct) + '\n')
            f.write("sbjct_start " + str(hsp.sbjct_start) + '\n')
    f.write("multiple_alignment " + str(blast_record.multiple_alignment) + '\n')
    f.close()


def online_alignment(seq, output=" "):
    try:
        print("Starting..............")
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)  # n becouse it is nucltide
        blast_record = NCBIXML.read(result_handle)
    except OSError as Error:
        print("Connect wrong ............... ")
        print(Error)
        return
    if output == " ":
        for description in blast_record.descriptions:
            print("******************************Description***************************************")
            print("Title ", description.title)
            print("Score ", description.score)
            print("e ", description.e)
            print("Number of Alignments ", description.num_alignments)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print("************************Alignments********************************")
                print("Title ", alignment.title)
                print("Length ", alignment.length)
                print("*****************************hsp********************************")
                print("score ", hsp.score)
                print("bits ", hsp.bits)
                print("expect ", hsp.expect)
                print("alignments ", hsp.num_alignments)
                print("identities ", hsp.identities)
                print("positives ", hsp.positives)
                print("gaps ", hsp.gaps)
                print("strand ", hsp.strand)
                print("frame ", hsp.frame)
                print("query ", hsp.query)
                print("query_start ", hsp.query_start)
                print("match ", hsp.match)
                print("sbjct ", hsp.sbjct)
                print("sbjct_start ", hsp.sbjct_start)
        print("multiple_alignment ", blast_record.multiple_alignment)
    else:
        write_blast(blast_record, output)
        print("the output of this alignment in the txt file----> ", output)
    print("Done.")


def merge_fasta(args):
    args = args.split(',')
    destination = input("enter the destination name of output file: ")
    m = open(destination, "w")
    meta = ''
    sequence = ''
    for arg in args:
        f = open(arg, 'r')
        for line in f:
            line = line.rstrip()
            if ">" in line:
                meta = line
                m.write("{} \n".format(meta))
                meta = ''
            else:
                sequence = sequence + line
                m.write("{} \n".format(sequence))
                sequence = ''
    print("Doneeeeeeeeee!!!")


def convert_to_fasta(gbk_filename):
    if os.path.exists(gbk_filename):
        if not (gbk_filename.endswith(".gb")):
            print("Please give the name of the file with the .gb extension!")
        else:
            faa_filename = input("enter fasta file name: ")
            if not (faa_filename.endswith(".fasta") or faa_filename.endswith(".fa")):
                print("Please give the name of the file with the .fasta extension!")
            else:
                with open(gbk_filename) as input_handle, open(faa_filename, "w") as output_handle:
                    sequences = SeqIO.parse(input_handle, "genbank")
                    count = SeqIO.write(sequences, output_handle, "fasta")
                print("Doneeeeeeeee!!!")
    else:
        print("invalid file !!")


def get_opt():
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "f:p:o:")
        for opt in range(len(opts)):
            if opts[opt][0] in ['-f']:
                if opts[opt][1] == 'gc_percentage':
                    if len(opts) == 2:  # opts = (-f , -p)
                        if opts[opt + 1][0] in ['-p']:
                            gc_percentage(opts[opt + 1][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("gc_percentage takes only one parameter")

                elif opts[opt][1] == 'dna_transcription':
                    if len(opts) == 2:
                        if opts[opt + 1][0] in ['-p']:
                            dna_transcription(opts[opt + 1][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("dna_transcription takes only one parameter")

                elif opts[opt][1] == 'rev_complement':
                    if len(opts) == 2:
                        if opts[opt + 1][0] in ['-p']:
                            rev_complement(opts[opt + 1][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("rev_complement takes only one parameter")

                elif opts[opt][1] == 'calc_nbases':
                    if len(opts) == 2:
                        if opts[opt + 1][0] in ['-p']:
                            calc_nbases(opts[opt + 1][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("calc_nbases takes only one parameter")

                elif opts[opt][1] == 'filter_nbases':
                    if len(opts) == 2:
                        if opts[opt + 1][0] in ['-p']:
                            filter_nbases(opts[opt + 1][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("filter_nbases takes only one parameter")

                elif opts[opt][1] == 'convert_to_fasta':
                    if len(opts) == 2:
                        if opts[opt + 1][0] in ['-p']:
                            convert_to_fasta(opts[opt + 1][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("convert_to_fasta takes one parameter")

                elif opts[opt][1] == 'is_valid':
                    if len(opts) == 3:
                        if opts[opt + 1][0] in ['-p']:
                            is_valid(opts[opt + 1][1], opts[opt + 2][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("is_valid takes two parameters")

                elif opts[opt][1] == 'seq_alignment':
                    if len(opts) == 4:
                        if opts[opt + 1][0] in ['-p'] and opts[opt + 3][0] in ['-o']:
                            seq_alignment(opts[opt + 1][1], opts[opt + 2][1], opts[opt + 3][1])
                        else:
                            print("wrong option!!!")
                    elif len(opts) == 3:
                        if opts[opt + 1][0] in ['-p']:
                            seq_alignment(opts[opt + 1][1], opts[opt + 2][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("Error in parameters!!")

                elif opts[opt][1] == 'seq_alignment_files':
                    if len(opts) == 4:
                        if opts[opt + 1][0] in ['-p'] and opts[opt + 3][0] in ['-o']:
                            seq_alignment_files(opts[opt + 1][1], opts[opt + 2][1], opts[opt + 3][1])
                        else:
                            print("wrong option!!!")
                    elif len(opts) == 3:
                        if opts[opt + 1][0] in ['-p']:
                            seq_alignment_files(opts[opt + 1][1], opts[opt + 2][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("Error in parameters!!")

                elif opts[opt][1] == 'online_alignment':
                    if len(opts) == 3:
                        if opts[opt + 1][0] in ['-p'] and opts[opt + 2][0] in ['-o']:
                            online_alignment(opts[opt + 1][1], opts[opt + 2][1])
                        else:
                            print("wrong option!!!")
                    elif len(opts) == 2:
                        if opts[opt + 1][0] in ['-p']:
                            online_alignment(opts[opt + 1][1])
                        else:
                            print("wrong option!!!")
                    else:
                        print("Error in parameters!!")

                elif opts[opt][1] == 'merge_fasta':
                    args = opts[opt + 1:]
                    arg_list = []
                    for arg in args:
                        arg_list.append(arg[1])
                    arg_string = ','.join(arg_list)
                    merge_fasta(arg_string)
    except getopt.GetoptError:
        print("Error, correct your input")


get_opt()

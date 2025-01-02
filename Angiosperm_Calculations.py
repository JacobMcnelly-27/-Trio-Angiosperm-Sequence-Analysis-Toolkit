#importing python libraries#
import os
import argparse
import sys


#define Main function of code including runtime arguments#
#run python3 Angiosperm_Calculations.py -h to see options#
def main():
    parser = argparse.ArgumentParser(description="A script that can compile and run simple statistical tests on aligned(you can align files using this script but you then must run any other function seperately)fasta file datasets of trios(annual,perrenial,outgroup)")
    parser.add_argument("-a","--annual", help ="This will be the fasta tag for the annual species of your trio")
    parser.add_argument("-p", "--perrenial", help="This will be the fasta tag for the perrenial species of your trio")
    parser.add_argument("-r", "--outgroup", help="This will be the fasta tag for the outgroup species of your trio")
    parser.add_argument("-t", "--type", help=" This will determine if you are running the tests on compilled or batch running sequences(c(compiled) or b(batch))")
    parser.add_argument("-o", "--output", help="This will be the assigned output directory or file depending on the function you are preforming")
    parser.add_argument("-i", "--input", help="This will be the assigned input directory or file depending on the function you are preforming")
    parser.add_argument("-f", "--function", help="This will be what function you wish to preform || align (aligns the directory of files chosen) || tajima (preforms Tajima's D1 test on the file or files chosen giving output as a CSV file with M1, M2, D1, D2, and chi-squared values) || phylip (Copies FASTA files and formats them as PHYLIP files to be used on server for PAML) ||")

    args = parser.parse_args()
    if args.function == "align":
        if args.type == "batch":
            for file in os.listdir(args.input):
                alignment(f"{args.input}/{file}",f"{args.output}/{file}")
        if args.type == "compiled":
            alignment(args.input, args.output)

    if args.function == "tajima":
        with open(f"{args.output}", "w") as f:
            f.write("filename,Da,Dp,M1,M2,chi^2\n")
        if args.type == "batch":
            for file in os.listdir(args.input):
                tajimas(args.output,args.input, file, args.annual, args.perrenial, args.outgroup)

        if args.type == "compiled":
            tajimas(args.output,"./",args.input, args.annual, args.perrenial, args.outgroup)

    if args.function == "phylip":
        if args.type == "batch":
            for file in os.listdir(args.input):
                PHYLIP(args.input, file, args.output)

        if args.type == "compiled":
            PHYLIP("./", args.input, "./")

def PHYLIP(input_dir,input_file,output):
    count = 0
    input = input_file.replace(".fasta", "")
    with open(f"{output}/{input}.phy", "w") as file:
        with open(f"{input_dir}/{input}.fasta", "r") as f:
            lines = f.readlines()
            num_seq = int(len(lines)/2)
            length_seq = len(lines[1].strip())
            file.write(f"{num_seq} {length_seq}\n")
            for line in lines:
                if ">" in line:
                    count += 1
                    line = line.replace(">", "")
                    line = line.replace("\n", "")
                    print(line)
                    length_name = len(line)
                    for x in range(0,10-length_name):
                        line += " "
                    line += lines[count]
                    file.write(f"{line}")
                else:
                    count += 1

def tajimas(output,input_dir, input_file, annual, perrenial, outgroup):
    annual_seq = aquire(annual, f"{input_dir}/{input_file}")
    perrenial_seq = aquire(perrenial, f"{input_dir}/{input_file}")
    outgroup_seq = aquire(outgroup, f"{input_dir}/{input_file}")

    gene = str(input_file)
    Dp = formula(differences(annual_seq, perrenial_seq), differences(perrenial_seq, outgroup_seq), differences(annual_seq, outgroup_seq))
    Da = formula(differences(annual_seq, perrenial_seq), differences(annual_seq, outgroup_seq), differences(perrenial_seq, outgroup_seq))
    M1 = m_vals(annual_seq, perrenial_seq, outgroup_seq)
    M2 = m_vals(perrenial_seq, annual_seq, outgroup_seq)
    chi = chi_square(M1, M2)
    print(f"{gene},{Dp},{Da},{M1},{M2},{chi}")
    with open(f"{output}", "a") as f:
        f.write(f"{gene},{Dp},{Da},{M1},{M2},{chi}\n")

#aquires sequence of specific species from provided aligned trio fasta file (either batch or compilled)#
def aquire(sID,file):
    seq_string = ""
    count = 0
    with open(f"./{file}", "r") as f:
        file_data = f.readlines()
        for line in file_data:
            if line == f">{sID}\n":
                count += 1
                try:
                    while '>' not in file_data[count]:
                        seq_string += file_data[count]
                        count += 1
                except IndexError:
                    pass
            count += 1
    return seq_string


#this funtion will count the differences of aligned sequences in a trio file(gaps will be skipped)#
def differences(seq_a,seq_b):
    diff = 0
    count = 0
    for char in seq_a:
        try:
            if char != seq_b[count]:
                if char != '-' and seq_b[count] != '-':
                    diff += 1
            else:
                pass
        except IndexError:
            pass
        count += 1
    return diff

#formula function simply allows for the usage of the formula for calculating Dp and Da from the Ohta 1993 EQ.5#
def formula(dpa,dpr,dar):
    #dp = (dpa + dpr - dar )/2#
    #da = (dpa + dar - dpr )/2#
    answer = (dpa + dpr - dar )/2
    return answer

#calcualte M1 and M2 values(m1 being number of differences only present in sequence 1 or the annual and m2 being thus for sequence 2 or the perrenial)#
#credit: Tajima 1993 Simple Methods for Testing the Molecular Evolutionary Clock Hypothesis#
def m_vals(seq_a,seq_b,seq_c):
    # m1 is when seq_a != seq_p and seq_a != seq_o
    # m2 is when seq_p != seq_a and seq_p != seq_o
    m_val = 0
    count = 0
    for char in seq_a:
        try:
            if char != seq_b[count] and char != seq_c[count]:
                if char != '-' and seq_b[count] != '-' and seq_c[count] != '-':
                    m_val += 1
                else:
                    pass
        except IndexError:
            pass
        count += 1
    return m_val

#calculate chi-squared values#
def chi_square(m1,m2):
    try:
        chi = (((m1 - m2)**2)/(m1 + m2))
    except ZeroDivisionError:
        chi = "negative number issue fix later!!!!!"
    return chi

##align sequences in order to compare##
def alignment(finalized_id,output):
    output.replace(".fasta", "")
    os.system(fr'java -jar ./macse_v2.07.jar -prog alignSequences -seq ./{finalized_id} -out_NT ./{output}_trans_align_NT.fasta -out_AA ./{output}_trans_align_AA.fasta')


if __name__ == "__main__":
    main()
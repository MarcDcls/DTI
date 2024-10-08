import os
import tqdm


##################################### PARAMETERS #####################################

# Target codons for the Deaminase with the index of the target nucleotide in the codon
# Example: ["CAG", 0] means that the target is the 1st nucleotide (C) of the codon CAG
targets = [["CAG", 0],
           ["CAA", 0],
           ["TGA", 1]]

# PAMs (Protospacer Adjacent Motifs)
# N represents any nucleotide
pams = ["NGG", "NGA"]

# Window for the nucleotide targetted by the deaminase before the PAM
deaminase_window = [15, 21]

# Threshold of concordance for a sequence to be considered as an off-target
concordance_threshold = 13

######################################################################################


def find_sequences(string: str, targets: list, pams: list, deaminase_window: list):
    """Display every possible sequence of nucleotides in the string beginning by 
    a target codon and having a PAM in the deaminase window.

    Args:
        string (str): sequence of nucleotides
        targets (list): list of target codons with the index of the target nucleotide in the codon
        pams (list): list of PAM sequences
        deaminase_window (list): upper and lower bounds of the deaminase window

    Returns:
        list: list of the possible sequences and their indexes in the string
    """ 
    results = []
    for i in range(len(string)//3):
        codon = string[i*3:i*3+3]
        for target in targets:
            if codon == target[0]:
                for j in range(deaminase_window[0]+target[1], deaminase_window[1]+1+target[1]):
                    if i*3+j+3 > len(string):
                        break
                    for pam in pams:
                        for k in range(3):
                            if pam[k] == "N":
                                continue
                            if string[i*3+j+k] != pam[k]:
                                break
                            if k == 2:
                                results.append({"index": i*3,
                                                "target": target,
                                                "pam": string[i*3+j:i*3+j+3],
                                                "sequence": string[i*3:i*3+j+3]})
    return results


def check_off_target(sequence: dict, genome: str, concordance_threshold: int):
    """Check if the sequence has off-targets in the genome sequence.

    Args:
        sequence (dict): sequence dict
        genome (str): genome sequence
        concordance_threshold (int): threshold of concordance for a sequence to be considered as an off-target
    
    Returns:
        list: list of off-targets index in the genome sequence
    """
    # Support for the circularity of the genome sequence
    genome = (genome + genome[:deaminase_window[1]+3]).upper()

    sequence_length = len(sequence["sequence"])
    off_targets = []
    for i in range(len(genome)-3):
        for k in range(3):
            if sequence["pam"][k] == "N":
                continue
            if genome[i+k] != sequence["pam"][k]:
                break
            if k == 2:
                concordance = 0
                for j in range(sequence_length-3):
                    if i-j-1 < 0:
                        break
                    if genome[i-j-1] == sequence["sequence"][sequence_length-j-4]:
                        concordance += 1
                    else:
                        break
                if concordance >= concordance_threshold:
                    off_targets.append({"index": i-sequence_length+3,
                                        "concordance": concordance,
                                        "sequence": genome[i-sequence_length+3:i+3]})
    return off_targets


def analyse_fasta_file(dir_path: str, fasta_file: str, genome_file: str, save: bool, targets: list, pams: list, deaminase_window: list, concordance_threshold: int):
    """Analyse a fasta file to find the possible sequences, excluding the off-targets.

    Args:
        dir_path (str): directory path
        fasta_file (str): fasta file name
        genome_file (str): genome file name
        save (bool): save the results in a file, if False the results are printed in the console
        targets (list): list of target codons with the index of the target nucleotide in the codon
        pams (list): list of PAM sequences
        deaminase_window (list): upper and lower bounds of the deaminase window
        concordance_threshold (int): threshold of concordance for a sequence to be considered as an off-target
    """
    # Open the output file if save is True
    if save:
        filename = fasta_file[:-6]
        output = open(dir_path + "/results_" + filename + ".txt", "w")
    else:
        output = None

    # Open the fasta file
    file_f = open(os.path.join(dir_path, fasta_file), "r")
    fasta = file_f.read()
    fasta = fasta.split("\n")

    # Open the genome file and format it
    file_g = open(os.path.join(dir_path, genome_file), "r")
    genome = file_g.read()
    genome = "".join(genome.split("\n")[1:])

    def write_output(text):
        if save:
            output.write(text + "\n")
        else:
            print(text)

    for i in tqdm.tqdm(range(len(fasta)//2)):
        write_output(fasta[i*2])
        write_output("")

        sample = fasta[i*2+1]
        sequences = find_sequences(sample, targets, pams, deaminase_window)
        for j, sequence in enumerate(sequences):
            write_output(f"Sequence {j+1}: " + sequence["sequence"])
            write_output("    PAM: " + sequence["pam"])
            write_output("    Target: " + sequence["target"][0])
            write_output("    Index in FASTA sample: " + str(sequence["index"]))
            write_output("")

            off_targets = check_off_target(sequence, genome, concordance_threshold)
            replicas = []
            for off_target in off_targets[:]:
                if off_target["sequence"] == sequence["sequence"]:
                    replicas.append(off_target)
                    off_targets.remove(off_target)

            write_output("    Number of exact replica: " + str(len(replicas)))
            for replica in replicas:
                write_output("        Index in genome: " + str(replica["index"]))
            write_output("")
            
            write_output("    Number of off-targets: " + str(len(off_targets)))
            write_output("")
            for off_target in off_targets:
                write_output("        Off-target sequence: " + off_target["sequence"])
                write_output("        Concordance: " + str(off_target["concordance"]))
                write_output("        Index in genome: " + str(off_target["index"]))
                write_output("")

        write_output("---------------------------------------------")
        write_output("")

    # Close the output file if save is True
    if save:
        output.close()


if __name__ == "__main__":

    for dir in os.listdir("data/"):
        dir_path = os.path.join("data", dir)
        
        if os.path.isdir(dir_path):
            results_files = [file for file in os.listdir(dir_path) if file.startswith("results_")]
            fasta_files = [file for file in os.listdir(dir_path) if file.endswith(".fasta")]
            genome_files = [file for file in os.listdir(dir_path) if file.endswith(".fa")]
            
            # Repository already analysed
            if len(results_files) == 1:
                continue

            # Check if the repository is correctly formatted
            if len(fasta_files) != 1 or len(genome_files) != 1:
                print(f"The repository {dir_path} should contain exactly one .fasta file and one .fa file !")
                break

            print(f"Analyzing {dir}...")
            analyse_fasta_file(dir_path, fasta_files[0], genome_files[0], True, targets, pams, deaminase_window, concordance_threshold)
            print(f"{dir} analysed !")
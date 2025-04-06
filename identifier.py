import os

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

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

# Number of missmatch allowed in an off-target
allowed_mismatches = 1

# Number of gaps allowed in an off-target
allowed_gaps = 1

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
                                                "pam": pam,
                                                "sequence": string[i*3:i*3+j+3]})
    return results

def step_off_target(off_targets, target, string, last_index, min_concordance, max_mismatches, max_gaps, concordance = 0, mismatches = 0, gaps = 0, result = ""):
    if len(string) == 0 or len(target) == 0:
        if concordance >= min_concordance and result[0] != "*" and result[0] != "-":
            off_targets.append({"index": [last_index-len(result), last_index],
                                "concordance": concordance,
                                "sequence": result,
                                "mismatches": mismatches,
                                "gaps": gaps})
        return
    
    if target[-1] == string[-1]:
        step_off_target(off_targets, target[:-1], string[:-1], last_index, min_concordance, max_mismatches, max_gaps, concordance + 1, mismatches, gaps, target[-1] + result)
    else:
        if mismatches < max_mismatches:
            step_off_target(off_targets, target[:-1], string[:-1], last_index, min_concordance, max_mismatches, max_gaps, concordance, mismatches + 1, gaps, "*" + result)
        if gaps < max_gaps:
            step_off_target(off_targets, target, string[:-1], last_index, min_concordance, max_mismatches, max_gaps, concordance, mismatches, gaps + 1, "-" + result)
        
        if concordance >= min_concordance and result[0] != "*" and result[0] != "-":
            off_targets.append({"index": [last_index-len(result), last_index],
                                "concordance": concordance,
                                "sequence": result,
                                "mismatches": mismatches,
                                "gaps": gaps})

def check_off_target(sequence: dict, genome: str, concordance_threshold: int, allowed_mismatches: int):
    """Check if the sequence has off-targets in the genome sequence.

    Args:
        sequence (dict): sequence dict
        genome (str): genome sequence
        concordance_threshold (int): threshold of concordance for a sequence to be considered as an off-target
    
    Returns:
        list: list of off-targets index in the genome sequence
    """
    # Support for the circularity of the genome sequence
    genome = (genome + genome[:2]).upper()

    sequence_length = len(sequence["sequence"])
    off_targets = []
    for i in range(len(genome)-2):
        for k in range(3):
            if sequence["pam"][k] == "N":
                continue
            if genome[i+k] != sequence["pam"][k]:
                break
            if k == 2:
                start = i - sequence_length + 3
                if start >= 0:
                    gen_seq = genome[start:i]
                else:
                    gen_seq = genome[start-2: -2] + genome[:i]
                step_off_target(off_targets, sequence["sequence"][:-3], gen_seq, i, concordance_threshold, allowed_mismatches, allowed_gaps, result=genome[i:i+3])

    return off_targets

def analyse_fasta_file(dir_path: str, fasta_file: str, genome_file: str, save: bool, targets: list, pams: list, deaminase_window: list, concordance_threshold: int, allowed_mismatches: int):
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
    buffer = ""

    # Open the fasta file
    file_f = open(os.path.join(dir_path, fasta_file), "r")
    fasta = file_f.read()
    fasta = fasta.split("\n")

    # Open the genome file and format it
    file_g = open(os.path.join(dir_path, genome_file), "r")
    genome = file_g.read()
    genome = "".join(genome.split("\n")[1:])
    CI_genome = get_CI_sequence(genome)
    size_genome = len(genome)

    # Analyse the fasta file
    nb_KO = 0
    nb_targets = 0
    for i in tqdm(range(len(fasta)//2)):
        buffer += fasta[i*2] + "\n\n"

        sample = fasta[i*2+1]
        sequences = find_sequences(sample, targets, pams, deaminase_window)

        buffer += f"Number of potential targets found: {len(sequences)}\n\n"
        nb_targets += len(sequences)

        if len(sequences) > 0:
            nb_KO += 1

        for j, sequence in enumerate(sequences):
            buffer += f"Target {j+1}: {sequence['sequence']}\n"
            buffer += f"    PAM: {sequence['pam']}\n"
            buffer += f"    Target: {sequence['target'][0]}\n"
            buffer += f"    Index in FASTA sample: {sequence['index']}\n\n"

            off_targets = check_off_target(sequence, genome, concordance_threshold, allowed_mismatches)
            replicas = []
            for off_target in off_targets[:]:
                if off_target["sequence"] == sequence["sequence"]:
                    replicas.append(off_target)
                    off_targets.remove(off_target)

            CI_off_targets = check_off_target(sequence, CI_genome, concordance_threshold, allowed_mismatches)
            CI_replicas = []
            for CI_off_target in CI_off_targets[:]:
                if CI_off_target["sequence"] == sequence["sequence"]:
                    CI_replicas.append(CI_off_target)
                    CI_off_targets.remove(CI_off_target)

            buffer += f"    Number of exact replica: {len(replicas) + len(CI_replicas)}\n\n"
            for replica in replicas:
                buffer += f"        Index in genome main strand: {replica['index'][0]} to {replica['index'][1]}\n"
            for CI_replica in CI_replicas:
                buffer += f"        Index in genome secondary strand: {size_genome - CI_replica['index'][0]} to {size_genome - CI_replica['index'][1]}\n"
            buffer += "\n"
            
            buffer += f"    Number of off-targets: {len(off_targets) + len(CI_off_targets)}\n\n"
            for off_target in off_targets:
                buffer += f"        Off-target sequence: {off_target['sequence']}\n"
                buffer += f"        Concordance: {off_target['concordance']}\n"
                buffer += f"        Index in genome main strand: {off_target['index'][0]} to {off_target['index'][1]}\n\n"
            for CI_off_target in CI_off_targets:
                buffer += f"        Off-target sequence: {CI_off_target['sequence']}\n"
                buffer += f"        Concordance: {CI_off_target['concordance']}\n"
                buffer += f"        Index in genome secondary strand: {size_genome - CI_off_target['index'][0]} to {size_genome - CI_off_target['index'][1]}\n\n"

        buffer += "---------------------------------------------\n\n"

    # Write statistics at the beginning of the file
    tmp = buffer
    buffer = f"Number of KO: {nb_KO}/{len(fasta)//2} - {round(nb_KO/(len(fasta)//2)*100, 2)}%\n"
    buffer += f"Number of potential targets found among all genes: {nb_targets}\n\n"
    buffer += "---------------------------------------------\n\n"
    buffer += tmp

    # Printing or saving the results
    if save:
        filename = fasta_file[:-6]
        output = open(dir_path + "/results_" + filename + ".txt", "w")
        output.write(buffer)
        output.close()
    else:
        print(buffer)


def get_CI_sequence(sequence: str):
    """Get the complementary inverse sequence of a DNA sequence.

    Args:
        sequence (str): DNA sequence

    Returns:
        str: complementary inverse DNA sequence
    """
    complementary = ""
    for nucleotide in sequence:
        if nucleotide == "A":
            complementary += "T"
        elif nucleotide == "a":
            complementary += "t"
        elif nucleotide == "T":
            complementary += "A"
        elif nucleotide == "t":
            complementary += "a"
        elif nucleotide == "C":
            complementary += "G"
        elif nucleotide == "c":
            complementary += "g"
        elif nucleotide == "G":
            complementary += "C"
        elif nucleotide == "g":
            complementary += "c"
        else:
            complementary += "N"
    return complementary[::-1]

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
            analyse_fasta_file(dir_path, fasta_files[0], genome_files[0], True, targets, pams, deaminase_window, concordance_threshold, allowed_mismatches)
            print(f"{dir} analysed !")
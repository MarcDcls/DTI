sequence_starts = ["caa", "cag"]
sequence_ends = ["gg"]
sequence_lengths = [15, 21]

# Display every sequence of codons (the sequence, its length and its beginning 
# index) of the gene placed at gene_filename having a beginning within starts, an 
# end within ends and a length in the range of lengths.
# If a result_filename is specified, save the result in this file.
def find_sequences(gene_filename, starts, ends, lengths, result_filename=None):
    file = open(gene_filename, "r")
    gene = file.read()
    result = ""

    sequence_count = 0
    for start in starts:
        sequence_start_index = 0
        while True:
            sequence_start_index = gene.find(start, sequence_start_index)
            if sequence_start_index == -1:
                break
            elif sequence_start_index%3 == 0:
                for end in ends:
                    sequence_end_index = gene.find(end, sequence_start_index+lengths[0]-2, sequence_start_index+lengths[1]+1)
                    if sequence_end_index != -1:
                        result += "\n--- SEQUENCE " + str(sequence_count+1) + " ---"
                        result += "\nstarting position : " + str(sequence_start_index)
                        result += "\nlength : " + str(sequence_end_index-sequence_start_index+2)
                        result += "\nsequence : " + gene[sequence_start_index:sequence_end_index+2]
                        result += "\n------------------\n"
                        sequence_count += 1
            sequence_start_index += 1

    if sequence_count == 0:
        result += "NO SEQUENCE FOUND"

    if result_filename != None:
        file = open(result_filename, "w")
        file.write(result)

    print(result)

find_sequences("gene.txt", sequence_starts, sequence_ends, sequence_lengths, "result.txt")
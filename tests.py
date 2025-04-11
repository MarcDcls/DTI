from identifier import *

if __name__ == "__main__":

    coding_target = ["CAG", 0, "coding"]
    coding_offset_target = ["CAG", 1, "coding"]
    complementary_target = ["TGA", 1, "complementary"]

    pams = ["NGG", "NGA"]

    deaminase_window = [15, 21]
    concordance_threshold = 13
    allowed_mismatches = 1
    allowed_gaps = 1
    nb_processes = cpu_count()

    ###############################################################

    find_test = True

    # Testing lower bound of the deaminase window
    seq = "ATGXXXXXXCAG123456789012NGG"
    if len(find_sequences(seq, coding_target[0], coding_target[1], coding_target[2], pams, deaminase_window)) != 1:
        print("Test 1 failed")
        find_test = False
    seq = "ATGXXXXXXCAG12345678901NGGX"
    if len(find_sequences(seq, coding_target[0], coding_target[1], coding_target[2], pams, deaminase_window)) != 0:
        print("Test 2 failed")
        find_test = False

    # Testing upper bound of the deaminase window
    seq = "ATGXXXXXXCAG123456789012345678NGG"
    if len(find_sequences(seq, coding_target[0], coding_target[1], coding_target[2], pams, deaminase_window)) != 1:
        print("Test 3 failed")
        find_test = False
    seq = "ATGXXXXXXCAG1234567890123456789NGGXX"
    if len(find_sequences(seq, coding_target[0], coding_target[1], coding_target[2], pams, deaminase_window)) != 0:
        print("Test 4 failed")
        find_test = False

    # Testing false positive on a shifted target
    seq = "ATGXXXXXCAG123456789012345NGG"
    if len(find_sequences(seq, coding_target[0], coding_target[1], coding_target[2], pams, deaminase_window)) != 0:
        print("Test 5 failed")
        find_test = False

    # Testing several targets
    seq = "ATGXXXXXXCAG123456789012NGGXXXCAG123456789012NGG"
    if len(find_sequences(seq, coding_target[0], coding_target[1], coding_target[2], pams, deaminase_window)) != 2:
        print("Test 6 failed")
        find_test = False
    
    # Testing every target/PAM combination
    for target in targets:
        for pam in pams:
            seq = "ATGXXXXXX" + target[0] + "123456789012345" + pam
            if len(find_sequences(seq, target[0], target[1], "coding", pams, deaminase_window)) != 1:
                print("Test 7 failed")
                find_test = False
            if find_sequences(seq, target[0], target[1], "coding", pams, deaminase_window)[0]["index"] != 9:
                print("Test 8 failed")
                find_test = False

    # Testing target nucleotide not in first position
    seq = "ATGXXXXXXCAG123456789012NGG"
    if len(find_sequences(seq, coding_offset_target[0], coding_offset_target[1], coding_offset_target[2], pams, deaminase_window)) != 0:
        print("Test 9 failed")
        find_test = False
    seq = "ATGXXXXXXCAG1234567890123NGGXX"
    if len(find_sequences(seq, coding_offset_target[0], coding_offset_target[1], coding_offset_target[2], pams, deaminase_window)) != 1:
        print("Test 10 failed")
        find_test = False
    seq = "ATGXXXXXXCAG1234567890123456789NGGXX"
    if len(find_sequences(seq, coding_offset_target[0], coding_offset_target[1], coding_offset_target[2], pams, deaminase_window)) != 1:
        print("Test 11 failed")
        find_test = False
    seq = "ATGXXXXXXCAG12345678901234567890NGGX"
    if len(find_sequences(seq, coding_offset_target[0], coding_offset_target[1], coding_offset_target[2], pams, deaminase_window)) != 0:
        print("Test 12 failed")
        find_test = False
    
    if find_test:
        print("All find_sequences tests passed")

    ###############################################################

    offtarget_test = True
    test_candidate = {"index": 9, "target": ["CAG", 0], "pam": "NGG", "sequence": None}

    # Testing circularity of the gene sequence
    test_candidate["sequence"] = "CAG123456789012345NGG"
    genome_test = "345NGGXXXXXXCAG123456789012345NGGXXXCAG123456789012"
    if len(check_off_target(test_candidate, genome_test, 13, 0, 0)) != 2:
        print("Test 13 failed")
        offtarget_test = False
    genome_test = "GGXXXXXXCAG123456789012345NGGXXXCAG123456789012345N"
    if len(check_off_target(test_candidate, genome_test, 13, 0, 0)) != 2:
        print("Test 14 failed")
        offtarget_test = False

    # Testing concordance threshold
    test_candidate["sequence"] = "CAG123456789012345NGG"
    genome_test = "345NGGXXXXXX456789012345NGGXXXCAG123456789012"
    if len(check_off_target(test_candidate, genome_test, 13, 0, 0)) != 1:
        print("Test 15 failed")
        offtarget_test = False
    genome_test = "345NGGXXXXXXXX3456789012345NGGXXXCAG123456789012"
    if len(check_off_target(test_candidate, genome_test, 13, 0, 0)) != 2:
        print("Test 16 failed")
        offtarget_test = False

    # Testing allowed mismatches
    genome_test = "345NGGXXXXXXXXXXXXXXXXXXNGGXXXCAG12345*789012"
    if len(check_off_target(test_candidate, genome_test, 13, 1, 0)) != 1:
        print("Test 17 failed")
        offtarget_test = False
    genome_test = "345NGGXXXXXXXXXXXXXXXXXXNGGXXXCAG12345*78*012"
    if len(check_off_target(test_candidate, genome_test, 13, 1, 0)) != 0:
        print("Test 18 failed")
        offtarget_test = False

    # Testing allowed gaps
    genome_test = "345NGGXXXXXXXXXXXXXXXXXXNGGXXXCAG12345-6789012"
    if len(check_off_target(test_candidate, genome_test, 13, 0, 1)) != 1:
        print("Test 19 failed")
        offtarget_test = False
    genome_test = "345NGGXXXXXXXXXXXXXXXXXXNGGXXXCAG12345-678-9012"
    if len(check_off_target(test_candidate, genome_test, 13, 0, 1)) != 0:
        print("Test 20 failed")
        offtarget_test = False
    
    if offtarget_test:
        print("All check_off_target tests passed")

    ###############################################################
    
    complementary_search_test = True

    CI_target = get_CI_sequence(complementary_target[0])

    seq = "ATGXXCCN1234567890123TGAXXX"
    if len(find_sequences(get_CI_sequence(seq), CI_target, 2 - complementary_target[1], complementary_target[2], pams, deaminase_window)) != 1:
        print("Test 21 failed")
        complementary_search_test = False

    seq = "ATGXXCCN123456789012TGAXXX"
    if len(find_sequences(get_CI_sequence(seq), CI_target, 2 - complementary_target[1], complementary_target[2], pams, deaminase_window)) != 0:
        print("Test 22 failed")
        complementary_search_test = False

    if complementary_search_test:
        print("All find_sequences complementary tests passed")    


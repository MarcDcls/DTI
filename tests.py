from identifier import *

if __name__ == "__main__":
    find_test = True

    # Testing lower bound of the deaminase window
    seq = "ATGXXXXXXCAG123456789012NGG"
    if len(find_sequences(seq, targets, pams, deaminase_window)) != 1:
        print("Test 1 failed")
        find_test = False
    seq = "ATGXXXXXXCAG12345678901NGGX"
    if len(find_sequences(seq, targets, pams, deaminase_window)) != 0:
        print("Test 2 failed")
        find_test = False

    # Testing upper bound of the deaminase window
    seq = "ATGXXXXXXCAG123456789012345678NGG"
    if len(find_sequences(seq, targets, pams, deaminase_window)) != 1:
        print("Test 3 failed")
        find_test = False
    seq = "ATGXXXXXXCAG1234567890123456789NGGXX"
    if len(find_sequences(seq, targets, pams, deaminase_window)) != 0:
        print("Test 4 failed")
        find_test = False

    # Testing false positive on a shifted target
    seq = "ATGXXXXXCAG123456789012345NGG"
    if len(find_sequences(seq, targets, pams, deaminase_window)) != 0:
        print("Test 5 failed")
        find_test = False

    # Testing several targets
    seq = "ATGXXXXXXCAG123456789012NGGXXXCAG123456789012NGG"
    if len(find_sequences(seq, targets, pams, deaminase_window)) != 2:
        print("Test 6 failed")
        find_test = False
    
    # Testing every target/PAM combination
    for target in targets:
        for pam in pams:
            seq = "ATGXXXXXX" + target[0] + "123456789012345" + pam
            if len(find_sequences(seq, targets, pams, deaminase_window)) != 1:
                print("Test 7 failed")
                find_test = False
            if find_sequences(seq, targets, pams, deaminase_window)[0]["index"] != 9:
                print("Test 8 failed")
                find_test = False

    # Testing target nucleotide not in first position
    test_targets = [("CAG", 1)]
    seq = "ATGXXXXXXCAG123456789012NGG"
    if len(find_sequences(seq, test_targets, pams, deaminase_window)) != 0:
        print("Test 9 failed")
        find_test = False
    seq = "ATGXXXXXXCAG1234567890123NGGXX"
    if len(find_sequences(seq, test_targets, pams, deaminase_window)) != 1:
        print("Test 10 failed")
        find_test = False
    seq = "ATGXXXXXXCAG1234567890123456789NGGXX"
    if len(find_sequences(seq, test_targets, pams, deaminase_window)) != 1:
        print("Test 11 failed")
        find_test = False
    seq = "ATGXXXXXXCAG12345678901234567890NGGX"
    if len(find_sequences(seq, test_targets, pams, deaminase_window)) != 0:
        print("Test 12 failed")
        find_test = False
    
    if find_test:
        print("All find_sequences tests passed")


    offtarget_test = True
    test_candidate = {"index": 9, "target": ["CAG", 0], "pam": "NGG", "sequence": None}

    # Testing circularity of the gene sequence
    test_candidate["sequence"] = "CAG123456789012345NGG"
    genome_test = "345NGGXXXXXXCAG123456789012345NGGXXXCAG123456789012"
    if len(check_off_target(test_candidate, genome_test, 13, 0)) != 2:
        print("Test 13 failed")
        offtarget_test = False
    genome_test = "GGXXXXXXCAG123456789012345NGGXXXCAG123456789012345N"
    if len(check_off_target(test_candidate, genome_test, 13, 0)) != 2:
        print("Test 14 failed")
        offtarget_test = False

    # Testing concordance threshold
    test_candidate["sequence"] = "CAG123456789012345NGG"
    genome_test = "345NGGXXXXXX456789012345NGGXXXCAG123456789012"
    if len(check_off_target(test_candidate, genome_test, 13, 0)) != 1:
        print("Test 15 failed")
        offtarget_test = False
    genome_test = "345NGGXXXXXXXX3456789012345NGGXXXCAG123456789012"
    if len(check_off_target(test_candidate, genome_test, 13, 0)) != 2:
        print("Test 16 failed")
        offtarget_test = False

    if offtarget_test:
        print("All check_off_target tests passed")
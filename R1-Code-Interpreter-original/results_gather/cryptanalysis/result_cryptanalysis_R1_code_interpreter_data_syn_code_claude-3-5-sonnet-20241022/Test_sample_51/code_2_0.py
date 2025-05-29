def check_all_conditions(n1, n2, l1, l2):
    def check_position_feedback(guess_n1, guess_n2, guess_l1, guess_l2, 
                              n_correct_pos, n_correct_wrong_pos, 
                              n_too_small, n_too_large,
                              l_correct_pos, l_correct_wrong_pos, l_too_early):
        # Check numbers
        numbers_correct_pos = sum([
            (guess_n1 == n1),
            (guess_n2 == n2)
        ])
        numbers_correct_wrong_pos = sum([
            (guess_n1 == n2),
            (guess_n2 == n1)
        ])
        numbers_too_small = sum([
            int(guess_n1) < min(int(n1), int(n2)),
            int(guess_n2) < min(int(n1), int(n2))
        ])
        numbers_too_large = sum([
            int(guess_n1) > max(int(n1), int(n2)),
            int(guess_n2) > max(int(n1), int(n2))
        ])
        
        # Check letters
        letters_correct_pos = sum([
            (guess_l1 == l1),
            (guess_l2 == l2)
        ])
        letters_correct_wrong_pos = sum([
            (guess_l1 == l2),
            (guess_l2 == l1)
        ])
        letters_too_early = sum([
            ord(guess_l1) < ord(l1) or ord(guess_l1) < ord(l2),
            ord(guess_l2) < ord(l1) or ord(guess_l2) < ord(l2)
        ])

        return (numbers_correct_pos == n_correct_pos and 
                numbers_correct_wrong_pos == n_correct_wrong_pos and
                numbers_too_small == n_too_small and
                numbers_too_large == n_too_large and
                letters_correct_pos == l_correct_pos and
                letters_correct_wrong_pos == l_correct_wrong_pos and
                (l_too_early == -1 or letters_too_early >= l_too_early))

    # Check all 10 conditions
    c1 = check_position_feedback('5', '2', 'F', 'R', 
                               0, 1, 1, 0, 0, 0, 2)
    c2 = check_position_feedback('0', '2', 'I', 'T',
                               0, 0, 2, 0, 0, 0, 2)
    c3 = check_position_feedback('3', '7', 'A', 'O',
                               0, 0, 0, 0, 0, 0, 2)
    c4 = check_position_feedback('1', '7', 'F', 'J',
                               0, 0, 0, 0, 0, 0, 2)
    c5 = check_position_feedback('7', '9', 'N', 'P',
                               0, 0, 0, 2, 0, 1, 1)
    c6 = check_position_feedback('9', '1', 'C', 'K',
                               0, 0, 0, 0, 0, 0, 2)
    c7 = check_position_feedback('5', '6', 'S', 'M',
                               0, 2, 0, 0, 0, 0, 2)
    c8 = check_position_feedback('7', '6', 'B', 'J',
                               0, 1, 0, 1, 0, 0, 2)
    c9 = check_position_feedback('3', '1', 'R', 'G',
                               0, 0, 2, 0, 0, 0, 2)
    c10 = check_position_feedback('9', '4', 'A', 'V',
                                0, 0, 0, 0, 1, 0, 1)

    return all([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10])

# Test all possible combinations with numbers 5 and 6
solutions = []
for n1 in ['5', '6']:
    for n2 in ['5', '6']:
        if n1 == n2:
            continue
        for l1 in ['N', 'P', 'V']:  # Limited to the letters we know are involved
            for l2 in ['N', 'P', 'V']:
                if l1 == l2:
                    continue
                if check_all_conditions(n1, n2, l1, l2):
                    solutions.append([n1, n2, l1, l2])

print(solutions)
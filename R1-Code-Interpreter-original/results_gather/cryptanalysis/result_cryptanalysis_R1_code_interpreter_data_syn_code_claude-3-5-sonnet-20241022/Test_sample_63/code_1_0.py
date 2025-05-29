def verify_guess(test_combo, guess_str, num_feedback, letter_feedback):
    # Split test combination and guess into components
    test_nums = [test_combo[0], test_combo[1]]
    test_letters = [test_combo[2], test_combo[3]]
    guess_nums = [int(guess_str[0]), int(guess_str[1])]
    guess_letters = [guess_str[2], guess_str[3]]
    
    # Check numbers
    correct_pos_n = 0
    wrong_pos_n = 0
    too_large_n = 0
    too_small_n = 0
    
    # First digit
    if test_nums[0] == guess_nums[0]:
        correct_pos_n += 1
    elif test_nums[0] == guess_nums[1]:
        wrong_pos_n += 1
    elif test_nums[0] > guess_nums[0]:
        too_large_n += 1
    else:
        too_small_n += 1
    
    # Second digit
    if test_nums[1] == guess_nums[1]:
        correct_pos_n += 1
    elif test_nums[1] == guess_nums[0]:
        wrong_pos_n += 1
    elif test_nums[1] > guess_nums[1]:
        too_large_n += 1
    else:
        too_small_n += 1
    
    # Check letters
    correct_pos_l = 0
    wrong_pos_l = 0
    too_early_l = 0
    too_late_l = 0
    
    # First letter
    if test_letters[0] == guess_letters[0]:
        correct_pos_l += 1
    elif test_letters[0] == guess_letters[1]:
        wrong_pos_l += 1
    elif ord(test_letters[0]) < ord(guess_letters[0]):
        too_early_l += 1
    else:
        too_late_l += 1
    
    # Second letter
    if test_letters[1] == guess_letters[1]:
        correct_pos_l += 1
    elif test_letters[1] == guess_letters[0]:
        wrong_pos_l += 1
    elif ord(test_letters[1]) < ord(guess_letters[1]):
        too_early_l += 1
    else:
        too_late_l += 1
    
    return [correct_pos_n, wrong_pos_n, too_large_n, too_small_n] == list(num_feedback) and \
           [correct_pos_l, wrong_pos_l, too_early_l, too_late_l] == list(letter_feedback)

# All test cases
test_cases = [
    ('38UW', (0,1,1,0), (0,0,0,2)),  # one number wrong pos, one too large, both letters too late
    ('69OG', (1,0,1,0), (1,0,1,0)),  # one number correct pos, one too large, one letter correct pos, one too early
    ('83OP', (1,0,1,0), (1,0,0,0)),  # one number correct pos, one too large, one letter correct pos
    ('48SP', (0,0,0,0), (0,0,0,0)),  # all incorrect
    ('42ZE', (0,0,0,0), (0,0,0,0)),  # all incorrect
    ('01WG', (0,0,0,2), (0,0,0,0)),  # both numbers too small
    ('07JX', (0,0,0,0), (0,0,0,0)),  # all incorrect
    ('31HK', (0,1,0,1), (0,0,2,0)),  # one number wrong pos, one too small, both letters too early
    ('37VY', (0,1,1,0), (0,0,0,2)),  # one number wrong pos, one too large, both letters too late
    ('62QH', (1,0,0,1), (0,0,0,0)),  # one number correct pos, one too small
    ('94CL', (0,0,0,0), (0,0,2,0)),  # both letters too early
    ('71TE', (0,0,0,0), (0,1,1,0))   # one letter wrong pos, one too early
]

# Try all possible combinations
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2:
                    continue
                combo = [n1, n2, l1, l2]
                valid = True
                for case in test_cases:
                    if not verify_guess(combo, case[0], case[1], case[2]):
                        valid = False
                        break
                if valid:
                    print(f"Found valid solution: {combo}")
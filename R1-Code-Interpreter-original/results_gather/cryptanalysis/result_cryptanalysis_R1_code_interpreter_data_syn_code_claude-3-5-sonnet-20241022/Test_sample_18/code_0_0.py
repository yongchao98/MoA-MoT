def check_guess(guess, actual):
    # Check if a guess matches all the given feedback conditions
    correct_pos_num = 0
    correct_wrong_pos_num = 0
    correct_pos_let = 0
    correct_wrong_pos_let = 0
    
    # Check numbers (first two positions)
    for i in range(2):
        if guess[i] == actual[i]:
            correct_pos_num += 1
        elif actual[i] in guess[:2]:
            correct_wrong_pos_num += 1
            
    # Check letters (last two positions)
    for i in range(2,4):
        if guess[i] == actual[i]:
            correct_pos_let += 1
        elif actual[i] in guess[2:]:
            correct_wrong_pos_let += 1
            
    return correct_pos_num, correct_wrong_pos_num, correct_pos_let, correct_wrong_pos_let

def is_later_in_alphabet(a, b):
    return a > b

def check_all_conditions(candidate):
    # Check against all given guesses and their feedback
    
    # Condition 1: 56QS
    c1_num_pos, c1_num_wrong, c1_let_pos, c1_let_wrong = check_guess("56QS", candidate)
    if not (c1_num_pos == 1 and c1_num_wrong == 0 and c1_let_pos == 0 and c1_let_wrong == 0):
        return False
    if not (is_later_in_alphabet('Q', candidate[2]) and is_later_in_alphabet('S', candidate[3])):
        return False
        
    # Condition 2: 47KB
    c2_num_pos, c2_num_wrong, c2_let_pos, c2_let_wrong = check_guess("47KB", candidate)
    if not (c2_num_pos == 0 and c2_num_wrong == 0 and c2_let_pos == 0 and c2_let_wrong == 0):
        return False
        
    # Condition 3: 83CN
    c3_num_pos, c3_num_wrong, c3_let_pos, c3_let_wrong = check_guess("83CN", candidate)
    if not (c3_num_pos == 0 and c3_num_wrong == 1 and c3_let_pos == 0 and c3_let_wrong == 1):
        return False
    if not int(candidate[0]) > 3 or not int(candidate[1]) > 3:  # at least one number should be > 3
        return False
    if not is_later_in_alphabet('N', candidate[3]):  # N is too late
        return False
        
    # Condition 4: 35JX
    c4_num_pos, c4_num_wrong, c4_let_pos, c4_let_wrong = check_guess("35JX", candidate)
    if not (c4_num_pos == 0 and c4_num_wrong == 1 and c4_let_pos == 0 and c4_let_wrong == 0):
        return False
    if not is_later_in_alphabet('J', candidate[2]) and not is_later_in_alphabet('X', candidate[3]):
        return False
        
    # Condition 5: 95FG
    c5_num_pos, c5_num_wrong, c5_let_pos, c5_let_wrong = check_guess("95FG", candidate)
    if not (c5_num_pos == 0 and c5_num_wrong == 1 and c5_let_pos == 0 and c5_let_wrong == 1):
        return False
        
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        candidate = n1 + n2 + l1 + l2
                        if check_all_conditions(candidate):
                            print([n1, n2, l1, l2])
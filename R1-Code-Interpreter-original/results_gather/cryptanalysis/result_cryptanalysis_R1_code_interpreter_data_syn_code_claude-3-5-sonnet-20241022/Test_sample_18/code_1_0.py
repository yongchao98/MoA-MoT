def check_guess(guess, actual):
    # Same as before
    correct_pos_num = 0
    correct_wrong_pos_num = 0
    correct_pos_let = 0
    correct_wrong_pos_let = 0
    
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            correct_pos_num += 1
        elif actual[i] in guess[:2]:
            correct_wrong_pos_num += 1
            
    # Check letters
    for i in range(2,4):
        if guess[i] == actual[i]:
            correct_pos_let += 1
        elif actual[i] in guess[2:]:
            correct_wrong_pos_let += 1
            
    return correct_pos_num, correct_wrong_pos_num, correct_pos_let, correct_wrong_pos_let

def check_all_conditions(candidate):
    # Condition 1: 56QS
    c1 = check_guess("56QS", candidate)
    if c1 != (1, 0, 0, 0):  # one number correct position, rest incorrect
        return False
    if not (ord(candidate[2]) < ord('Q') and ord(candidate[3]) < ord('S')):
        return False
    
    # Condition 2: 47KB
    c2 = check_guess("47KB", candidate)
    if c2 != (0, 0, 0, 0):  # all incorrect
        return False
    
    # Condition 3: 83CN
    c3 = check_guess("83CN", candidate)
    if c3 != (0, 1, 0, 1):  # one number and one letter wrong position
        return False
    if not (int(candidate[0]) > 3 and int(candidate[1]) > 3):  # numbers should be > 3
        return False
    if not ord(candidate[3]) < ord('N'):  # letter should be before N
        return False
    
    # Condition 4: 35JX
    c4 = check_guess("35JX", candidate)
    if c4 != (0, 1, 0, 0):  # one number wrong position, both letters incorrect
        return False
    if not (ord(candidate[2]) < ord('J') and ord(candidate[3]) < ord('X')):
        return False
    
    # Condition 5: 95FG
    c5 = check_guess("95FG", candidate)
    if c5 != (0, 1, 0, 1):  # one number and one letter wrong position
        return False
    
    # Additional constraints from analyzing the conditions:
    # From condition 1: 5 must be in correct position (first position)
    if candidate[0] != '5':
        return False
    
    # From condition 3: 8 must be one of the numbers (but wrong position)
    if '8' not in candidate[:2]:
        return False
    
    # From condition 5: F or G must be one of the letters (but wrong position)
    if not ('F' in candidate[2:] or 'G' in candidate[2:]):
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
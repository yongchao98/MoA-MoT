def check_guess(candidate):
    n1, n2, l1, l2 = candidate
    
    # Helper function to check if a letter is in wrong position
    def letter_in_wrong_pos(letter, pos1, pos2, guess1, guess2):
        return (letter == guess1 and pos1 != guess2) or (letter == guess2 and pos2 != guess1)
    
    # Condition 1: 04XW - both numbers wrong, letters too late
    if n1 == '0' or n2 == '4' or l1 >= 'X' or l2 >= 'W':
        return False
    
    # Condition 2: 06IH - numbers wrong, one letter correct but wrong pos, one too late
    if n1 == '0' or n2 == '6':
        return False
    has_h_wrong_pos = letter_in_wrong_pos('H', l1, l2, 'I', 'H')
    if not has_h_wrong_pos:
        return False
    
    # Condition 3: 54FA - all incorrect
    if n1 == '5' or n2 == '4' or l1 == 'F' or l2 == 'A':
        return False
    
    # Condition 4: 07ES - numbers wrong, one letter correct but wrong pos
    if n1 == '0' or n2 == '7':
        return False
    has_e_wrong_pos = letter_in_wrong_pos('E', l1, l2, 'E', 'S')
    if not has_e_wrong_pos:
        return False
    
    # Condition 5: 08HI - one number correct in pos, one too small
    correct_num_pos = (n1 == '0' and n2 != '8') or (n1 != '0' and n2 == '8')
    if not correct_num_pos:
        return False
    
    # Condition 6: 57OR - numbers wrong, letters too late
    if n1 == '5' or n2 == '7' or l1 >= 'O' or l2 >= 'R':
        return False
    
    # Condition 7: 29FB - one number correct in pos, one too large
    correct_num_pos_2 = (n1 == '2' and n2 != '9') or (n1 != '2' and n2 == '9')
    if not correct_num_pos_2:
        return False
    if l1 == 'F' or l2 == 'B':
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        if check_guess([n1, n2, l1, l2]):
                            valid_solutions.append([n1, n2, l1, l2])

print(valid_solutions)
# Possible numbers and letters based on deductions
possible_numbers = [1, 8]
possible_letters = ['O', 'X']

# Verify the combination
combination = [8, 1, 'O', 'X']

# Check against each condition
def check_combination(combination):
    # Condition 1: 95VB
    if (combination[0] in [9, 5] and combination[0] != 9) or (combination[1] in [9, 5] and combination[1] != 5):
        if combination[2] not in ['V', 'B'] and combination[3] not in ['V', 'B']:
            # Condition 2: 34OH
            if combination[0] not in [3, 4] and combination[1] not in [3, 4]:
                if combination[2] == 'O' and combination[3] not in ['H']:
                    # Condition 3: 63TK
                    if combination[0] not in [6, 3] and combination[1] not in [6, 3]:
                        if combination[2] not in ['T', 'K'] and combination[3] not in ['T', 'K']:
                            # Condition 4: 61XN
                            if (combination[0] in [6, 1] and combination[0] != 6) or (combination[1] in [6, 1] and combination[1] != 1):
                                if combination[2] == 'X' and combination[3] not in ['N']:
                                    return True
    return False

# Check if the combination is correct
is_correct = check_combination(combination)
print(is_correct)
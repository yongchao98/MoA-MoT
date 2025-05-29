def check_guess(letters):
    # Condition 1: 28MV
    # - one letter correct and in correct position
    # - one letter incorrect and too early in alphabet
    c1_correct_count = (letters[0] == 'M') + (letters[1] == 'V')
    if c1_correct_count != 1:  # Exactly one letter must be correct
        return False
    
    # For the other letter (not the correct one), it must be too early
    if letters[0] == 'M':  # If M is correct
        if not (letters[1] < 'V'):  # Other letter must be before V
            return False
    elif letters[1] == 'V':  # If V is correct
        if not (letters[0] < 'M'):  # Other letter must be before M
            return False
    else:
        return False

    # Condition 2: 41FX
    # Both letters must be different
    if 'F' in letters or 'X' in letters:
        return False

    # Condition 3: 58VP
    # - one letter correct but wrong position
    # - one letter incorrect and too early
    c3_wrong_pos = ('V' in letters and letters[1] != 'V') or ('P' in letters and letters[0] != 'P')
    if not c3_wrong_pos:
        return False
    # The other letter must be too early
    if 'V' in letters:
        if not (letters[1] < 'P'):
            return False
    elif 'P' in letters:
        if not (letters[0] < 'V'):
            return False
    else:
        return False

    # Condition 4: 93IT
    # - one letter correct but wrong position
    # - one letter incorrect and too early
    c4_wrong_pos = ('I' in letters and letters[1] != 'I') or ('T' in letters and letters[0] != 'T')
    if not c4_wrong_pos:
        return False
    # The other letter must be too early
    if 'I' in letters:
        if not (letters[1] < 'T'):
            return False
    elif 'T' in letters:
        if not (letters[0] < 'I'):
            return False
    else:
        return False

    return True

# Test all possible combinations
solutions = []
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
for l1 in letters:
    for l2 in letters:
        if l1 != l2:
            if check_guess([l1, l2]):
                solutions.append([2, 8, l1, l2])

print("Solutions:", solutions)
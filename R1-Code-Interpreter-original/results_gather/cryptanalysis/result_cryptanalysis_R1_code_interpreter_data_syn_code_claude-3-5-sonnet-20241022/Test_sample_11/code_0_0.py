import string

def is_too_early(letter1, letter2):
    return letter1 < letter2

def check_guess(numbers, letters):
    # Check against all conditions
    
    # Condition 1: 28MV
    correct_pos_count = 0
    if letters[0] == 'M':
        correct_pos_count += 1
    if letters[1] == 'V':
        correct_pos_count += 1
    if correct_pos_count != 1:
        return False
    # One letter should be too early
    early_count = 0
    for l in letters:
        if l < 'M' or l < 'V':
            early_count += 1
    if early_count != 1:
        return False

    # Condition 2: 41FX
    if 'F' in letters or 'X' in letters:
        return False

    # Condition 3: 58VP
    v_pos = letters.index('V') if 'V' in letters else -1
    p_pos = letters.index('P') if 'P' in letters else -1
    if v_pos == -1 and p_pos == -1:
        return False
    if v_pos != -1 and p_pos != -1:
        return False
    if v_pos == 1 or p_pos == 1:  # Wrong position check
        return False
    early_count = 0
    for l in letters:
        if l < 'V' or l < 'P':
            early_count += 1
    if early_count != 1:
        return False

    # Condition 4: 93IT
    i_pos = letters.index('I') if 'I' in letters else -1
    t_pos = letters.index('T') if 'T' in letters else -1
    if i_pos == -1 and t_pos == -1:
        return False
    if i_pos != -1 and t_pos != -1:
        return False
    if i_pos == 1 or t_pos == 1:  # Wrong position check
        return False
    early_count = 0
    for l in letters:
        if l < 'I' or l < 'T':
            early_count += 1
    if early_count != 1:
        return False

    return True

# Generate all possible letter combinations
possible_letters = []
for l1 in string.ascii_uppercase:
    for l2 in string.ascii_uppercase:
        if l1 != l2:
            if check_guess([2, 8], [l1, l2]):
                possible_letters.append([l1, l2])

print([2, 8] + possible_letters[0])
def check_guess(letters):
    # Helper function to check if a letter is too early compared to another
    def is_too_early(l, compared_to):
        return l < compared_to and l != compared_to
    
    # Condition 1: 28MV
    # - one letter correct and in correct position
    # - one letter incorrect and too early
    c1_correct_pos = (letters[0] == 'M') != (letters[1] == 'V')  # exactly one correct
    c1_early = sum(1 for l in letters if (l != 'M' and is_too_early(l, 'M')) or 
                                       (l != 'V' and is_too_early(l, 'V'))) == 1
    if not (c1_correct_pos and c1_early):
        return False

    # Condition 2: 41FX
    # - both letters must be different from F and X
    if 'F' in letters or 'X' in letters:
        return False

    # Condition 3: 58VP
    # - one letter correct but wrong position
    # - one letter incorrect and too early
    c3_wrong_pos = sum(1 for l in ['V', 'P'] if l in letters) == 1  # exactly one present
    if 'V' in letters and letters[1] == 'V':  # V in wrong position
        return False
    if 'P' in letters and letters[0] == 'P':  # P in wrong position
        return False
    c3_early = sum(1 for l in letters if l not in ['V', 'P'] and 
                   (is_too_early(l, 'V') or is_too_early(l, 'P'))) == 1
    if not (c3_wrong_pos and c3_early):
        return False

    # Condition 4: 93IT
    # - one letter correct but wrong position
    # - one letter incorrect and too early
    c4_wrong_pos = sum(1 for l in ['I', 'T'] if l in letters) == 1  # exactly one present
    if 'I' in letters and letters[1] == 'I':  # I in wrong position
        return False
    if 'T' in letters and letters[0] == 'T':  # T in wrong position
        return False
    c4_early = sum(1 for l in letters if l not in ['I', 'T'] and 
                   (is_too_early(l, 'I') or is_too_early(l, 'T'))) == 1
    if not (c4_wrong_pos and c4_early):
        return False

    return True

# Generate and test all possible letter combinations
solutions = []
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
for l1 in letters:
    for l2 in letters:
        if l1 != l2:
            if check_guess([l1, l2]):
                solutions.append([2, 8, l1, l2])

print("Solutions:", solutions)
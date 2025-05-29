def check_guess(letters):
    # Helper function to count matching letters in correct positions
    def count_matches(guess_letters, test_letters):
        return sum(1 for i in range(2) if guess_letters[i] == test_letters[i])

    # Helper function to count letters that are too early
    def count_too_early(letter, reference_letters):
        return any(letter < ref and letter != ref for ref in reference_letters)

    # Condition 1: 28MV
    # - Both numbers correct (already confirmed as 2,8)
    # - One letter correct and in correct position
    # - One letter incorrect and too early
    matches_mv = count_matches(letters, ['M', 'V'])
    if matches_mv != 1:  # Exactly one letter must match position
        return False
    
    # Count letters that are too early compared to M or V
    early_letters_mv = sum(1 for l in letters if l not in ['M', 'V'] and (l < 'M' or l < 'V'))
    if early_letters_mv != 1:  # Exactly one letter must be too early
        return False

    # Condition 2: 41FX
    # - Both letters must be different from F and X
    if 'F' in letters or 'X' in letters:
        return False

    # Condition 3: 58VP
    # - One letter correct but in wrong position
    # - One letter incorrect and too early
    vp_in_wrong_pos = sum(1 for l in ['V', 'P'] if l in letters and letters.index(l) != ['V', 'P'].index(l))
    if vp_in_wrong_pos != 1:  # Exactly one letter must be present but in wrong position
        return False
    
    # Count letters that are too early compared to V or P
    non_vp_letters = [l for l in letters if l not in ['V', 'P']]
    early_letters_vp = sum(1 for l in non_vp_letters if l < 'V' or l < 'P')
    if early_letters_vp != 1:  # Exactly one non-V/P letter must be too early
        return False

    # Condition 4: 93IT
    # - One letter correct but in wrong position
    # - One letter incorrect and too early
    it_in_wrong_pos = sum(1 for l in ['I', 'T'] if l in letters and letters.index(l) != ['I', 'T'].index(l))
    if it_in_wrong_pos != 1:  # Exactly one letter must be present but in wrong position
        return False
    
    # Count letters that are too early compared to I or T
    non_it_letters = [l for l in letters if l not in ['I', 'T']]
    early_letters_it = sum(1 for l in non_it_letters if l < 'I' or l < 'T')
    if early_letters_it != 1:  # Exactly one non-I/T letter must be too early
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
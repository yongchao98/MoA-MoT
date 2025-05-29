def verify_guess(l1, l2):
    # Condition 1: 28MV
    # - One letter (M or V) must be in correct position
    # - The other letter must be too early in alphabet
    mv_correct = (l1 == 'M') != (l2 == 'V')  # exactly one correct
    if not mv_correct:
        return False
    
    # Check if non-matching letter is too early
    if l1 == 'M':  # M is correct in first position
        if not (l2 < 'V'):  # second letter should be before V
            return False
    elif l2 == 'V':  # V is correct in second position
        if not (l1 < 'M'):  # first letter should be before M
            return False

    # Condition 2: 41FX
    # - Neither F nor X can be in the password
    if 'F' in (l1, l2) or 'X' in (l1, l2):
        return False

    # Condition 3: 58VP
    # - One letter (V or P) must be in wrong position
    # - Other letter must be too early
    vp_wrong_pos = (l1 == 'V' and 'V' in (l1, l2)) or (l2 == 'P' and 'P' in (l1, l2))
    if not vp_wrong_pos:
        return False
    
    # Check if the other letter is too early
    other_letter = l2 if l1 in 'VP' else l1
    if not (other_letter < 'P'):
        return False

    # Condition 4: 93IT
    # - One letter (I or T) must be in wrong position
    # - Other letter must be too early
    it_wrong_pos = (l1 == 'I' and 'I' in (l1, l2)) or (l2 == 'T' and 'T' in (l1, l2))
    if not it_wrong_pos:
        return False
    
    # Check if the other letter is too early
    other_letter = l2 if l1 in 'IT' else l1
    if not (other_letter < 'I'):
        return False

    return True

# Test all possible combinations
solutions = []
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
for l1 in letters:
    for l2 in letters:
        if l1 != l2:
            if verify_guess(l1, l2):
                solutions.append([2, 8, l1, l2])

print("All solutions:", solutions)
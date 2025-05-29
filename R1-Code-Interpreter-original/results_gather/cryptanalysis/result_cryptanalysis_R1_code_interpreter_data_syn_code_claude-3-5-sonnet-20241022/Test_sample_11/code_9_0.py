def verify_guess(l1, l2):
    # Condition 1: 28MV
    # - One letter must be correct and in correct position
    # - The other letter must be too early in alphabet
    mv_correct = (l1 == 'M') + (l2 == 'V')
    if mv_correct != 1:  # Exactly one must be correct
        return False
    
    # Check if the non-matching letter is too early
    if l1 == 'M':  # If M is correct
        if not (l2 < 'V'):  # Second letter must be before V
            return False
    elif l2 == 'V':  # If V is correct
        if not (l1 < 'M'):  # First letter must be before M
            return False

    # Condition 2: 41FX
    # Both letters must be different from F and X
    if 'F' in (l1, l2) or 'X' in (l1, l2):
        return False

    # Condition 3: 58VP
    # - One letter must be present but in wrong position
    # - The other letter must be too early
    vp_present = ('V' in (l1, l2)) + ('P' in (l1, l2))
    if vp_present != 1:  # Exactly one must be present
        return False
    
    # Check wrong position
    if 'V' in (l1, l2) and l2 == 'V':  # V can't be in second position
        return False
    if 'P' in (l1, l2) and l1 == 'P':  # P can't be in first position
        return False
    
    # Check if other letter is too early
    if 'V' in (l1, l2):
        other = l2 if l1 == 'V' else l1
        if not (other < 'P'):
            return False
    else:  # P is present
        other = l1 if l2 == 'P' else l2
        if not (other < 'V'):
            return False

    # Condition 4: 93IT
    # - One letter must be present but in wrong position
    # - The other letter must be too early
    it_present = ('I' in (l1, l2)) + ('T' in (l1, l2))
    if it_present != 1:  # Exactly one must be present
        return False
    
    # Check wrong position
    if 'I' in (l1, l2) and l2 == 'I':  # I can't be in second position
        return False
    if 'T' in (l1, l2) and l1 == 'T':  # T can't be in first position
        return False
    
    # Check if other letter is too early
    if 'I' in (l1, l2):
        other = l2 if l1 == 'I' else l1
        if not (other < 'T'):
            return False
    else:  # T is present
        other = l1 if l2 == 'T' else l2
        if not (other < 'I'):
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
def verify_guess(l1, l2):
    # Condition 1: 28MV
    if not ((l1 == 'M' and l2 != 'V') or (l2 == 'V' and l1 != 'M')):
        return False
    
    if l1 == 'M':  # M in first position
        if l2 >= 'V':  # second letter must be before V
            return False
    else:  # V in second position
        if l1 >= 'M':  # first letter must be before M
            return False

    # Condition 2: 41FX
    if 'F' in (l1, l2) or 'X' in (l1, l2):
        return False

    # Condition 3: 58VP
    vp_count = ('V' in (l1, l2)) + ('P' in (l1, l2))
    if vp_count != 1:
        return False
    
    if 'V' in (l1, l2):
        if l2 == 'V':  # V must be in first position
            return False
        if l2 >= 'P':  # other letter must be before P
            return False
    else:  # P is present
        if l1 == 'P':  # P must be in second position
            return False
        if l1 >= 'V':  # other letter must be before V
            return False

    # Condition 4: 93IT
    it_count = ('I' in (l1, l2)) + ('T' in (l1, l2))
    if it_count != 1:
        return False
    
    if 'I' in (l1, l2):
        if l2 == 'I':  # I must be in first position
            return False
        if l2 >= 'T':  # other letter must be before T
            return False
    else:  # T is present
        if l1 == 'T':  # T must be in second position
            return False
        if l1 >= 'I':  # other letter must be before I
            return False

    return True

solutions = []
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
for l1 in letters:
    for l2 in letters:
        if l1 != l2:
            if verify_guess(l1, l2):
                solutions.append([2, 8, l1, l2])

print("All solutions:", solutions)
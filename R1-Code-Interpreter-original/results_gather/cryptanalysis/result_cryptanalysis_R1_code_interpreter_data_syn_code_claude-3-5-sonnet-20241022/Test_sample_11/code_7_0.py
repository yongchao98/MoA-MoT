def verify_guess(l1, l2):
    # Condition 1: 28MV
    # "one letter is correct and in the correct position"
    c1_correct_pos = (l1 == 'M' or l2 == 'V')  # one must be correct
    if not c1_correct_pos:
        return False
    if l1 == 'M' and l2 == 'V':  # but not both
        return False
    
    # "one letter is incorrect and too early in the alphabet"
    if l1 == 'M':  # if M is correct
        if not (l2 < 'V'):  # other letter must be before V
            return False
    elif l2 == 'V':  # if V is correct
        if not (l1 < 'M'):  # other letter must be before M
            return False

    # Condition 2: 41FX
    # "both letters are incorrect"
    if 'F' in (l1, l2) or 'X' in (l1, l2):
        return False

    # Condition 3: 58VP
    # "one letter is correct but in the wrong position"
    c3_wrong_pos = ('V' in (l1, l2) and l2 != 'V') or ('P' in (l1, l2) and l1 != 'P')
    if not c3_wrong_pos:
        return False
    
    # "one letter is incorrect and too early"
    if 'V' in (l1, l2):
        other = l2 if l1 == 'V' else l1
        if not (other < 'P'):
            return False
    elif 'P' in (l1, l2):
        other = l1 if l2 == 'P' else l2
        if not (other < 'V'):
            return False

    # Condition 4: 93IT
    # "one letter is correct but in the wrong position"
    c4_wrong_pos = ('I' in (l1, l2) and l2 != 'I') or ('T' in (l1, l2) and l1 != 'T')
    if not c4_wrong_pos:
        return False
    
    # "one letter is incorrect and too early"
    if 'I' in (l1, l2):
        other = l2 if l1 == 'I' else l1
        if not (other < 'T'):
            return False
    elif 'T' in (l1, l2):
        other = l1 if l2 == 'T' else l2
        if not (other < 'I'):
            return False

    return True

# Test all possible combinations
solutions = []
for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
        if l1 != l2 and verify_guess(l1, l2):
            solutions.append([2, 8, l1, l2])

print("All solutions:", solutions)
def check_conditions(n1, n2, l1, l2):
    # Basic validations
    if sorted([n1, n2]) != sorted(['5', '6']):  # From condition 7
        return False
    
    if 'V' not in [l1, l2]:  # From condition 10
        return False
        
    if not ('N' in [l1, l2] or 'P' in [l1, l2]):  # From condition 5
        return False

    # Condition 1: 52FR
    # - one number correct but wrong position (5)
    # - one number too small (2)
    if '5' in [n1, n2]:
        if [n1, n2].index('5') == 0:  # if 5 is in first position
            return False

    # Condition 5: 79NP
    # - both numbers too large
    # - one letter correct but wrong position
    if 'N' in [l1, l2]:
        if [l1, l2].index('N') == 1:  # N should be in wrong position
            return False

    # Condition 8: 76BJ
    # - one number correct but wrong position (6)
    # - one number too large (7)
    if '6' in [n1, n2]:
        if [n1, n2].index('6') == 1:  # if 6 is in second position
            return False

    # Condition 10: 94AV
    # - both numbers wrong
    # - one letter correct and in correct position (V)
    if 'V' in [l1, l2]:
        if [l1, l2].index('V') != 0:  # V must be in first position
            return False

    return True

# Test all possible combinations
solutions = []
for n1, n2 in [('5','6'), ('6','5')]:
    for l1 in ['V']:  # V must be in first position
        for l2 in ['N', 'P']:  # N or P must be in second position
            if check_conditions(n1, n2, l1, l2):
                solutions.append([n1, n2, l1, l2])

print("Solutions found:", solutions)
def verify_combination(cand):
    # Clue 1: 54VU
    # - both numbers too small (✓ for all)
    # - one letter correct but wrong position
    # - one letter incorrect
    guess1_letters = set(['V', 'U'])
    cand_letters = set(cand[2:])
    if len(guess1_letters & cand_letters) != 1:  # must have exactly one match
        return False
    
    # Clue 2: 19LA
    # - both numbers incorrect (✓ for all)
    # - one letter correct but wrong position
    # - one letter incorrect
    guess2_letters = set(['L', 'A'])
    if len(guess2_letters & cand_letters) != 1:  # must have exactly one match
        return False
    
    # Clue 3: 38JV
    # - one number correct but wrong position (8 is correct but wrong pos)
    # - one number too small (3 is too small)
    # - one letter correct and in correct position
    # - one letter incorrect
    if cand[0] != '8':  # 8 must be first number
        return False
    if 'V' != cand[2]:  # V must be in first letter position
        return False
    
    # Clue 4: 14EC
    # - both numbers too small (✓ for all)
    # - both letters incorrect (✓ for all)
    
    # Clue 5: 47QU
    # - one number correct and in correct position (7)
    # - one number too small (4)
    # - both letters incorrect
    if cand[1] != '7':  # 7 must be second number
        return False
    
    return True

candidates = [[8, 7, 'V', 'A'], [8, 7, 'V', 'L'], [8, 7, 'L', 'V'], [8, 7, 'A', 'V']]
for cand in candidates:
    if verify_combination(cand):
        print(cand)
def check_letter_position(l1, l2):
    # Clue 1: KJ - both incorrect
    if l1 in 'KJ' or l2 in 'KJ':
        return False
    
    # Clue 2: SC - one correct wrong pos, one incorrect
    sc_count = sum(1 for l in [l1, l2] if l in 'SC')
    if sc_count != 1:  # exactly one letter should be either S or C
        return False
    
    # Clue 3: OQ - both incorrect
    if l1 in 'OQ' or l2 in 'OQ':
        return False
    
    # Clue 4: BV - one correct in position
    if l1 != 'B':  # B must be in first letter position
        return False
    if l2 == 'V':  # V is too late in alphabet
        return False
    
    return True

# Try all possible letter combinations
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
solution = None

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # letters can't repeat
            if check_letter_position(l1, l2):
                solution = ['3', '6', l1, l2]
                
print(solution)
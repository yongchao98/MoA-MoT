import string

def check_conditions(l1, l2):
    # From guess 1 (18VC)
    # One letter correct and in position, one too early
    correct_pos = (l1 == 'V' or l2 == 'C')
    if not correct_pos:
        return False
    
    # From guess 2 (54OF)
    # Both letters too early means our letters should be significantly after O,F
    if not (l1 > 'T' or l2 > 'T'):  # Using T as threshold (several letters after O)
        return False
    
    # From guess 3 (36AV)
    # V must be present but in wrong position from AV
    if 'V' not in [l1, l2]:
        return False
    
    # From guess 4 (42GE)
    # Both letters too early means our letters should be significantly after G,E
    if not (l1 > 'T' or l2 > 'T'):  # Using T as threshold
        return False
    
    # From guess 5 (08NQ)
    if l1 in 'NQ' or l2 in 'NQ':
        return False
    
    # Must be different letters
    if l1 == l2:
        return False
    
    # When a letter is "too early", the correct letter should be
    # significantly later in the alphabet (at least 5 letters after)
    if l1 == 'V':  # V is correct in first position
        if l2 <= 'T':  # Second letter should be well after the "too early" letters
            return False
    
    return True

valid_combinations = []

for l1 in string.ascii_uppercase:
    for l2 in string.ascii_uppercase:
        if check_conditions(l1, l2):
            valid_combinations.append(['3', '6', l1, l2])

print(valid_combinations)
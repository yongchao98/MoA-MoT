import string

def check_conditions(l1, l2):
    # From guess 1 (18VC)
    # One letter correct and in position, one too early
    correct_pos = (l1 == 'V' or l2 == 'C')
    if not correct_pos:
        return False
    
    # From guess 2 (54OF)
    # Both letters must be later than O and F respectively
    if l1 <= 'O' or l2 <= 'F':
        return False
    
    # From guess 3 (36AV)
    # V must be present but in wrong position from AV
    if 'V' not in [l1, l2]:
        return False
    if l1 <= 'A' or l2 <= 'A':
        return False
    
    # From guess 4 (42GE)
    # Both letters must be later than G and E
    if l1 <= 'G' or l2 <= 'E':
        return False
    
    # From guess 5 (08NQ)
    if l1 in 'NQ' or l2 in 'NQ':
        return False
    
    # Must be different letters
    if l1 == l2:
        return False
    
    # Additional constraint from guess 1:
    # If V is in position 1, the second letter must be later than C
    if l1 == 'V' and l2 <= 'C':
        return False
    # If V is in position 2, the first letter must be later than V
    if l2 == 'V' and l1 <= 'V':
        return False
    
    # From guess 4 (42GE): both letters too early
    # This means our letters must be LATER than both G and E
    if l1 <= 'G' or l2 <= 'G':
        return False
    
    return True

valid_combinations = []

for l1 in string.ascii_uppercase:
    for l2 in string.ascii_uppercase:
        if check_conditions(l1, l2):
            valid_combinations.append(['3', '6', l1, l2])

print(valid_combinations)
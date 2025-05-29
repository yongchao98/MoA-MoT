import string

def check_conditions(l1, l2):
    # From guess 1 (18VC)
    # One letter correct and in position, one too early
    correct_pos_count = 0
    if l1 == 'V' or l2 == 'C':  # Check if either letter is in correct position
        correct_pos_count += 1
    if correct_pos_count != 1:  # Must have exactly one correct position
        return False
    
    # The other letter must be later in alphabet than its corresponding position
    if l1 == 'V':  # If V is correct, C's position must have later letter
        if not l2 > 'C':
            return False
    elif l2 == 'C':  # If C is correct, V's position must have later letter
        if not l1 > 'V':
            return False
    
    # From guess 2 (54OF)
    # Both letters must be later than O and F respectively
    if not (l1 > 'O' and l2 > 'F'):
        return False
    
    # From guess 3 (36AV)
    # V must be present but in wrong position, other letter must be later than A
    has_v = (l1 == 'V' or l2 == 'V')
    if not has_v:
        return False
    if l1 == 'A' or l2 == 'A':  # A is too early
        return False
    
    # From guess 4 (42GE)
    # Both letters must be later than G and E respectively
    if not (l1 > 'G' and l2 > 'E'):
        return False
    
    # From guess 5 (08NQ)
    # Neither letter can be N or Q
    if l1 in 'NQ' or l2 in 'NQ':
        return False
    
    # Letters must be different
    if l1 == l2:
        return False
    
    return True

valid_combinations = []
letters = string.ascii_uppercase

for l1 in letters:
    for l2 in letters:
        if check_conditions(l1, l2):
            valid_combinations.append(['3', '6', l1, l2])

print(valid_combinations)
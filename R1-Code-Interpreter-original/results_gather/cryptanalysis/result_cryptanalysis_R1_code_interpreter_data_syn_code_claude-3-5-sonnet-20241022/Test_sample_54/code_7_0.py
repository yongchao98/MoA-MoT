def check_guess(candidate):
    n1, n2, l1, l2 = candidate
    
    # Numbers must be 2,8
    if n1 != '2' or n2 != '8':
        return False
    
    # Condition 1 (04XW): letters must be before W
    if l1 >= 'W' or l2 >= 'W':
        return False
    
    # Condition 2 (06IH): one letter must be H in wrong position
    if 'H' not in [l1, l2]:
        return False
    
    # Condition 4 (07ES): E must be present
    if 'E' not in [l1, l2]:
        return False
    
    # Condition 5 (08HI): one letter must be correct and in position
    if l1 != 'H' and l2 != 'I':
        return False
    
    # Condition 6 (57OR): letters before O
    if l1 >= 'O' or l2 >= 'O':
        return False
    
    # Condition 7 (29FB): letters can't be F or B
    if 'F' in [l1, l2] or 'B' in [l1, l2]:
        return False
    
    # Additional constraints:
    # - Letters can't be A (from condition 3)
    if 'A' in [l1, l2]:
        return False
    
    # - Letters must be different
    if l1 == l2:
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # letters can't repeat
            if check_guess(['2', '8', l1, l2]):
                valid_solutions.append(['2', '8', l1, l2])

print(valid_solutions)
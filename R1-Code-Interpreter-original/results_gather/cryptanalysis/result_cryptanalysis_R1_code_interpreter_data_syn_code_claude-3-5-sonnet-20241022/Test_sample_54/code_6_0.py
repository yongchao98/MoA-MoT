def check_guess(candidate):
    n1, n2, l1, l2 = candidate
    
    # Numbers must be 2,8 (already confirmed)
    if n1 != '2' or n2 != '8':
        return False
    
    # From condition 2 (06IH): H must be present but NOT in last position
    if 'H' in [l1, l2]:
        if l2 == 'H':  # H can't be in last position
            return False
    
    # From condition 4 (07ES): E must be present but NOT in last position
    if 'E' in [l1, l2]:
        if l2 == 'E':  # E can't be in last position
            return False
    
    # From condition 5 (08HI): Either H or I must be in same position as in 08HI
    if l1 != 'H' and l2 != 'I':  # Must have either H in third or I in fourth
        return False
    
    # Letters must be before O
    if l1 >= 'O' or l2 >= 'O':
        return False
    
    # Letters can't be F, A, or B
    if 'F' in [l1, l2] or 'A' in [l1, l2] or 'B' in [l1, l2]:
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
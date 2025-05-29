def check_position_constraints(guess):
    n1, n2, l1, l2 = guess
    
    # From guess 2 (10FX): if 0 is present, it must be in wrong position
    if 0 in [n1, n2]:
        if n2 == 0:  # 0 is in second position in guess 2
            return False
    
    # From guess 4 (48IM): I and M must be in wrong positions
    if l1 == 'I' and l2 == 'M':  # Same positions as guess 4
        return False
    
    # From guess 5 (86NM): M must be in wrong position
    if l2 == 'M':  # M in same position as guess 5
        return False
    
    return True

# Try all possible combinations
valid_solutions = []
for n1 in range(4):
    for n2 in range(4):
        if n1 == n2:
            continue
        for l1 in ['I', 'M']:
            for l2 in ['I', 'M']:
                if l1 == l2:
                    continue
                guess = [n1, n2, l1, l2]
                if check_guess(guess[:2], guess[2:], None) and check_position_constraints(guess):
                    valid_solutions.append(guess)
print(valid_solutions)
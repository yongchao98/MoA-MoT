def check_constraints(guess):
    n1, n2, l1, l2 = guess
    
    # Check number constraints
    if n1 > 5 or n2 > 5:  # Numbers can't be too large
        return False
    if n1 == n2:  # Numbers can't repeat
        return False
    
    # From clue 4: one number is 8 (but in wrong position)
    if n1 != 8 and n2 != 8:
        return False
        
    # Check letter constraints
    if l1 == l2:  # Letters can't repeat
        return False
    
    # Letters that are too late in alphabet
    late_letters = set('TKLZWQ')
    if l1 in late_letters or l2 in late_letters:
        return False
        
    # G must be in correct position (from clue 7)
    if l1 != 'G':
        return False
        
    # A must be present but in wrong position (from clue 6)
    if l2 != 'A':
        return False
    
    # E and O are incorrect
    if l1 in 'EO' or l2 in 'EO':
        return False
    
    return True

# Generate all possible combinations
valid_numbers = list(range(0, 6))  # 0-5
valid_letters = [c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' if c not in 'TKLZWQEO']

solutions = []
for n1 in valid_numbers:
    for n2 in valid_numbers:
        for l1 in valid_letters:
            for l2 in valid_letters:
                guess = [n1, n2, l1, l2]
                if check_constraints(guess):
                    solutions.append(guess)

print(solutions)
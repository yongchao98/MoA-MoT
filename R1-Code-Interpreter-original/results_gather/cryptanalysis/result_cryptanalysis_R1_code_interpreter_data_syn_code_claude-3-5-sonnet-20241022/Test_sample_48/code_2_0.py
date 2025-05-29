def is_valid_solution(n1, n2, l1, l2):
    numbers = [n1, n2]
    letters = [l1, l2, 'E']
    
    # Previous conditions remain the same...
    if 5 in numbers or 7 in numbers or 0 in numbers or 1 in numbers or 2 in numbers:
        return False
    if l1 in ['R', 'U', 'T', 'K', 'G', 'A'] or l2 in ['R', 'U', 'T', 'K', 'G', 'A']:
        return False
    
    # Additional specific conditions:
    
    # From guess 2 (26ZE):
    # - 6 must be in the solution (as it's the only valid number from this guess)
    if 6 not in numbers:
        return False
    
    # From guess 6 (94KA):
    # - Since 6 must be in solution, 9 must be the correct number (4 would conflict)
    if 9 not in numbers:
        return False
    
    # From guess 5 (25XY):
    # - Y must be in a different position than in this guess (where it was last)
    # - X must be in the solution as it's the only possible letter that works
    if l1 != 'X' and l2 != 'X':
        return False
    if 'Y' not in [l1, l2]:
        return False
    
    # The numbers must be [9,6] based on all conditions
    if sorted(numbers) != [6, 9]:
        return False
    
    return True

# Generate all possible combinations
solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in ['X', 'Y']:
            for l2 in ['X', 'Y']:
                if l1 != l2:
                    if is_valid_solution(n1, n2, l1, l2):
                        solutions.append([n1, n2, l1, l2, 'E'])

print(solutions)
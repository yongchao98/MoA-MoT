def is_valid_solution(n1, n2, l1, l2):
    numbers = [n1, n2]
    letters = [l1, l2, 'E']
    
    # Previous core conditions
    if sorted(numbers) != [6, 9]:  # Must be 6 and 9
        return False
    if set([l1, l2]) != set(['X', 'Y']):  # Must be X and Y
        return False
    
    # From guess 2 (26ZE):
    # - 6 must be in wrong position compared to this guess
    # - Since 6 is in position 2 here, it must be in position 1 in solution
    if n2 == 6:
        return False
    
    # From guess 6 (94KA):
    # - 9 must be in wrong position compared to this guess
    # - Since 9 is in position 1 here, it must be in position 2 in solution
    if n1 == 9:
        return False
    
    # From guess 5 (25XY):
    # - Y was in wrong position (position 4)
    # - Therefore Y cannot be in position 2 (which would be position 4 counting from 1)
    if l2 == 'Y':
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
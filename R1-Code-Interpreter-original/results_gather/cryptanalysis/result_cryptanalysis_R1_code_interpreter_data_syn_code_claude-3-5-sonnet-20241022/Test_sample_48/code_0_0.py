def check_guess(numbers, letters):
    # Check against all known conditions
    
    # Condition 1: 57RU
    if 5 in numbers or 7 in numbers:
        return False
    
    # Condition 2: 26ZE
    correct_count = sum(x in numbers for x in [2, 6])
    if correct_count != 1:
        return False
    if 'E' not in letters:
        return False
    if letters.index('E') != 3:  # E must be in last position
        return False
    
    # Condition 3: 01TK
    if 0 in numbers or 1 in numbers:
        return False
    
    # Condition 4: 56GR
    correct_count = sum(x in numbers for x in [5, 6])
    if correct_count != 1:
        return False
    
    # Condition 5: 25XY
    if 2 in numbers or 5 in numbers:
        return False
    if 'Y' not in letters:
        return False
    if letters.index('Y') == letters.index('E'):  # Y can't be in same position as E
        return False
    
    # Condition 6: 94KA
    correct_count = sum(x in numbers for x in [9, 4])
    if correct_count != 1:
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        numbers = [n1, n2]
        
        # Based on conditions, we know Y and E are in the solution
        # Try different positions for Y (except last position where E must be)
        for y_pos in range(2):  # Y can only be in first two positions
            letters = [None, None]
            letters[y_pos] = 'Y'
            letters.append('E')
            
            if check_guess(numbers, letters):
                valid_solutions.append([str(n1), str(n2), letters[0], letters[1], letters[2]])

print(valid_solutions)
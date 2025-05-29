def check_numbers(n1, n2):
    # Check against all conditions for numbers
    if n1 == n2:  # Numbers can't repeat
        return False
    
    # Condition 1: 25 both incorrect
    if n1 in [2, 5] or n2 in [2, 5]:
        return False
    
    # Condition 2: 28 both incorrect
    if n1 in [2, 8] or n2 in [2, 8]:
        return False
    
    # Condition 3: 93 - one correct in position, one too large
    correct_9_3 = (n1 == 9 and n2 != 3) or (n2 == 3 and n1 != 9)
    if not correct_9_3:
        return False
    
    # Condition 4: 24 both incorrect
    if n1 in [2, 4] or n2 in [2, 4]:
        return False
    
    # Condition 5: 07 - one correct but wrong position, one too small
    if not ((n1 == 7 and n2 != 0) or (n2 == 0 and n1 != 7)):
        return False
    
    return True

def check_letters(l1, l2):
    # Check against all conditions for letters
    if l1 == l2:  # Letters can't repeat
        return False
    
    # Known solution: DM in correct positions
    if l1 != 'D' or l2 != 'M':
        return False
    
    return True

# Find all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if check_numbers(n1, n2):
            if check_letters('D', 'M'):
                valid_solutions.append([str(n1), str(n2), 'D', 'M'])

print(valid_solutions)
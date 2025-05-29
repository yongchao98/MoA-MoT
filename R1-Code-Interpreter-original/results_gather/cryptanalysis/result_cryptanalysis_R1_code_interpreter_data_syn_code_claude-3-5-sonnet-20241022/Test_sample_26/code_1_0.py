def check_conditions(nums, letters):
    n1, n2 = nums
    l1, l2 = letters
    
    # Convert to strings for comparison
    solution = [str(n1), str(n2), l1, l2]
    
    # Condition 1: 83VM
    # Both numbers wrong, one letter correct but wrong position, one letter wrong
    if n1 == 8 or n2 == 3:
        return False
    if (l1 == 'V' and l2 == 'M') or (l1 == 'M' and l2 == 'V'):
        return False
    if not ((l1 == 'V' and l2 != 'M') or (l1 == 'M' and l2 != 'V') or 
            (l2 == 'V' and l1 != 'M') or (l2 == 'M' and l1 != 'V')):
        return False
    
    # Condition 2: 40WL
    # Both numbers and letters wrong
    if n1 == 4 or n2 == 0 or l1 in 'WL' or l2 in 'WL':
        return False
    
    # Condition 3: 50VD
    # Both numbers wrong, one letter correct but wrong position, one letter too early
    if n1 == 5 or n2 == 0:
        return False
    if 'D' in [l1, l2]:
        return False
    if not ('V' in [l1, l2]):
        return False
    
    # Condition 4: 40MC
    # Both numbers and letters wrong
    if n1 == 4 or n2 == 0 or l1 in 'MC' or l2 in 'MC':
        return False
    
    # Condition 5: 18JC
    # One number correct and in position, one number too large, both letters wrong
    correct_num = False
    if n1 == 1:
        correct_num = True
    elif n2 == 8:
        return False
    if not correct_num:
        return False
    if l1 in 'JC' or l2 in 'JC':
        return False
    
    # Condition 6: 41TW
    # One number correct but wrong position, one number wrong, both letters wrong
    if not (n1 == 4 or n2 == 1):
        return False
    if l1 in 'TW' or l2 in 'TW':
        return False
    
    # Condition 7: 07BF
    # Both numbers wrong, one letter correct but wrong position, one letter too early
    if n1 == 0 or n2 == 7:
        return False
    if 'B' in [l1, l2]:
        return False
    if 'F' in [l1, l2] and ('F' == l1):  # F is too early
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 != n2:
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 != l2:
                        if check_conditions((n1, n2), (l1, l2)):
                            valid_solutions.append([n1, n2, l1, l2])

print(valid_solutions)
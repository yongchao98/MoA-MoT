def is_valid_solution(n1, n2, l1, l2):
    numbers = [n1, n2]
    letters = [l1, l2, 'E']  # We know E is the last letter
    
    # Condition 1: 57RU - both numbers wrong, letters too early
    if 5 in numbers or 7 in numbers:
        return False
    if l1 in ['R', 'U'] or l2 in ['R', 'U']:
        return False
    
    # Condition 2: 26ZE - one number correct but wrong pos, one too small
    c2_nums = sum(x in [2, 6] for x in numbers)
    if c2_nums != 1:  # exactly one of 2 or 6 should be present
        return False
    if 2 in numbers:  # 2 is too small based on feedback
        return False
    
    # Condition 3: 01TK - both numbers too small, letters too early
    if 0 in numbers or 1 in numbers:
        return False
    if l1 in ['T', 'K'] or l2 in ['T', 'K']:
        return False
    
    # Condition 4: 56GR - one number correct but wrong pos, one too small
    c4_nums = sum(x in [5, 6] for x in numbers)
    if c4_nums != 1:
        return False
    if 5 in numbers:  # 5 is too small based on feedback
        return False
    
    # Condition 5: 25XY - both numbers too small, Y correct but wrong pos
    if 2 in numbers or 5 in numbers:
        return False
    if 'Y' not in [l1, l2]:  # Y must be present but not in last position
        return False
    
    # Condition 6: 94KA - one number correct but wrong pos, one too small
    c6_nums = sum(x in [9, 4] for x in numbers)
    if c6_nums != 1:
        return False
    if l1 in ['K', 'A'] or l2 in ['K', 'A']:
        return False
    
    return True

# Generate all possible combinations
solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        # Try all possible letters (after E in alphabet)
        for l1 in ['X', 'Y', 'Z']:
            for l2 in ['X', 'Y', 'Z']:
                if l1 != l2 and (l1 == 'Y' or l2 == 'Y'):  # Y must be present
                    if is_valid_solution(n1, n2, l1, l2):
                        solutions.append([n1, n2, l1, l2, 'E'])

print(solutions)
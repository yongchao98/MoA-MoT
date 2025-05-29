def check_combination(comb):
    n1, n2, l1, l2 = comb
    
    # Condition 1: "86AF" - both numbers incorrect, both letters incorrect
    if n1 in [8,6] or n2 in [8,6] or l1 in ['A','F'] or l2 in ['A','F']:
        return False
        
    # Condition 2: "98LF" - both numbers incorrect and too large
    if n1 in [9,8] or n2 in [9,8] or l1 in ['L','F'] or l2 in ['L','F']:
        return False
        
    # Condition 3: "20XK" - one number correct and in position, one too small
    correct_pos = (n1 == 2 and n2 != 0) or (n2 == 0 and n1 != 2)
    if not correct_pos or l1 in ['X','K'] or l2 in ['X','K']:
        return False
        
    # Condition 4: "38TM" - both numbers incorrect, both letters incorrect
    if n1 in [3,8] or n2 in [3,8] or l1 in ['T','M'] or l2 in ['T','M']:
        return False
        
    # Condition 5: "81JK" - both numbers incorrect, both letters incorrect
    if n1 in [8,1] or n2 in [8,1] or l1 in ['J','K'] or l2 in ['J','K']:
        return False
        
    # Condition 6: "23VQ" - one number correct and in position
    correct_pos = (n1 == 2 and n2 != 3) or (n2 == 3 and n1 != 2)
    if not correct_pos or l1 in ['V','Q'] or l2 in ['V','Q']:
        return False
        
    # Condition 7: "85IC" - both numbers incorrect, both letters incorrect
    if n1 in [8,5] or n2 in [8,5] or l1 in ['I','C'] or l2 in ['I','C']:
        return False
        
    # Condition 8: "03OE" - both numbers incorrect, both letters incorrect
    if n1 in [0,3] or n2 in [0,3] or l1 in ['O','E'] or l2 in ['O','E']:
        return False
        
    # Condition 9: "78WK" - one number correct but wrong position, one letter correct in position
    if not ((n1 == 8 and n2 == 7) or (n1 == 7 and n2 == 8)):
        return False
    if not (l1 == 'W' or l2 == 'K'):
        return False
        
    # Condition 10: "68KC" - both numbers incorrect, both letters incorrect
    if n1 in [6,8] or n2 in [6,8] or l1 in ['K','C'] or l2 in ['K','C']:
        return False
        
    # Condition 11: "47OM" - one number correct and in position
    correct_pos = (n1 == 4 and n2 != 7) or (n2 == 7 and n1 != 4)
    if not correct_pos or l1 in ['O','M'] or l2 in ['O','M']:
        return False
        
    # Condition 12: "58ZQ" - both numbers incorrect, both letters incorrect
    if n1 in [5,8] or n2 in [5,8] or l1 in ['Z','Q'] or l2 in ['Z','Q']:
        return False
        
    # Condition 13: "49DW" - both numbers incorrect, one letter correct but wrong position
    if n1 in [4,9] or n2 in [4,9]:
        return False
    if not ((l1 == 'W' and l2 != 'D') or (l2 == 'W' and l1 != 'D')):
        return False
        
    return True

# Test all possible combinations
solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2:
                    continue
                if check_combination((n1, n2, l1, l2)):
                    solutions.append([n1, n2, l1, l2])

print(f"Found {len(solutions)} solution(s):")
for sol in solutions:
    print(sol)
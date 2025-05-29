def verify_solution(sol):
    n1, n2, l1, l2 = sol
    
    # Condition 5: 18JC - confirms n1 must be 1
    if n1 != 1:
        return False
        
    # Condition 6: 41TW - confirms n2 must be 4
    if n2 != 4:
        return False
    
    # Condition 1: 83VM - V must be present but in different position than M
    if 'V' not in [l1, l2]:
        return False
        
    # Condition 3: 50VD - V must be present, D too early
    if 'D' in [l1, l2]:
        return False
        
    # Condition 7: 07BF - F must be in correct position (not too early)
    if 'F' in [l1, l2]:
        if l1 == 'F':  # F can't be first letter as it's "too early"
            return False
        
    return True

solutions = [[1, 2, 'F', 'V'], [1, 4, 'F', 'V'], [1, 5, 'F', 'V'], [1, 6, 'F', 'V'], [1, 9, 'F', 'V']]
valid = [s for s in solutions if verify_solution(s)]
print(valid)
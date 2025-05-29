def check_solution(n1, n2, l1, l2):
    # Condition 1: 52FR
    # - one number correct but wrong position
    # - one number too small
    # - both letters wrong
    c1 = ((n1 == '2' or n2 == '2') and 
          ('5' in [n1, n2]) and 
          l1 not in 'FR' and l2 not in 'FR')

    # Condition 2: 02IT
    # - both numbers too small
    c2 = int(n1) > 2 and int(n2) > 2

    # Condition 5: 79NP
    # - both numbers too large
    # - one letter correct but wrong position
    c5 = (int(n1) < 7 and int(n2) < 7 and
          ('N' in [l1, l2] or 'P' in [l1, l2]))

    # Condition 7: 56SM
    # - both numbers correct but wrong positions
    c7 = sorted([n1, n2]) == sorted(['5', '6'])

    # Condition 8: 76BJ
    # - one number correct but wrong position
    # - one number too large
    c8 = ('6' in [n1, n2] and '7' not in [n1, n2])

    # Condition 10: 94AV
    # - both numbers wrong
    # - one letter correct and in correct position
    c10 = ('V' in [l1, l2])

    return all([c1, c2, c5, c7, c8, c10])

# Test all possible combinations
solutions = []
for n1, n2 in [('5','6'), ('6','5')]:
    for l1 in ['N', 'P', 'V']:
        for l2 in ['N', 'P', 'V']:
            if l1 != l2:
                if check_solution(n1, n2, l1, l2):
                    solutions.append([n1, n2, l1, l2])

print("All possible solutions:", solutions)

# Additional verification
def verify_solution(sol):
    n1, n2, l1, l2 = sol
    
    # Verify V is in correct position (from condition 10)
    if 'V' not in [l1, l2]:
        return False
    
    # Verify N/P position (from condition 5)
    if 'N' in [l1, l2] and 'P' not in [l1, l2]:
        return True
    
    return False

final_solutions = [sol for sol in solutions if verify_solution(sol)]
print("\nFinal solution:", final_solutions[0] if final_solutions else "No solution found")
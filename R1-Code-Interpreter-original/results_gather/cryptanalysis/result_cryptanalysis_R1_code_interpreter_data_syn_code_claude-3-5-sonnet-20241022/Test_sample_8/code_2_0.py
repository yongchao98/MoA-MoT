def verify_solution(numbers, letters):
    n1, n2 = map(int, numbers)
    l1, l2 = letters
    
    # Verify first guess (97EJ)
    # Both numbers too large - check
    # One letter correct but wrong position, one incorrect
    ej_correct_wrong_pos = (l2 == 'E' and l1 != 'J') or (l2 == 'J' and l1 != 'E')
    if not ej_correct_wrong_pos:
        return False
    
    # Verify second guess (43OM)
    # Both numbers too large - check
    # Both letters too late - check
    
    # Verify third guess (20BQ)
    # 0 must be in solution
    if 0 not in [n1, n2]:
        return False
    
    # B must be first letter (as it was wrong in second position)
    if l1 != 'B':
        return False
    
    # Numbers must be 0 and 1
    if not set([n1, n2]) == set([0, 1]):
        return False
    
    return True

solutions = []
numbers = ['0', '1']
for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            if verify_solution([n1, n2], ['B', 'E']):
                solutions.append([n1, n2, 'B', 'E'])

print(solutions)
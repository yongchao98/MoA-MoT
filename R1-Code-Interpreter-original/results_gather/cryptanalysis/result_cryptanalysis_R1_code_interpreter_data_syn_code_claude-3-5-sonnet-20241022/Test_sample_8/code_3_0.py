def final_verify(solution):
    n1, n2, l1, l2 = solution
    
    # First guess (97EJ)
    if not (l2 == 'E' and l1 == 'B'):  # E must be second letter, B first letter
        return False
    
    # Third guess (20BQ)
    if not (n1 == '0'):  # 0 must be first number (as it was wrong in second position)
        return False
    
    if not (n2 == '1'):  # 1 is the only valid second number
        return False
    
    return True

solutions = [['0', '1', 'B', 'E'], ['1', '0', 'B', 'E']]
final_solution = [sol for sol in solutions if final_verify(sol)]
print(final_solution)
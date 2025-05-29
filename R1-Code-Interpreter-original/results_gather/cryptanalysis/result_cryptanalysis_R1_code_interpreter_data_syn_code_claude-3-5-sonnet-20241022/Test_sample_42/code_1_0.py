def check_letter_positions(solution):
    # Get positions of X and W in the solution
    x_pos = solution.index('X')
    w_pos = solution.index('W')
    
    # From guess 5 (98XJ): X cannot be in position 2
    if x_pos == 2:
        return False
        
    # From guess 6 (37XW): 
    # If X is in position 3, it's in same position as guess 6 (wrong)
    # If W is in position 3, it's in same position as guess 6 (wrong)
    if x_pos == 3 or w_pos == 3:
        return False
    
    return True

solutions = [[1, 2, 'W', 'X'], [1, 2, 'X', 'W'], [6, 2, 'W', 'X'], [6, 2, 'X', 'W']]
valid = []
for sol in solutions:
    if check_letter_positions(sol):
        valid.append(sol)
        
print(valid)
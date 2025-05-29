def final_check(solution):
    # From guess 4 (28JC):
    # 2 is correct but wrong position, meaning it can't be where it was in guess 4 (position 0)
    if solution[0] == 2:  # 2 can't be in first position
        return False
        
    # From guess 5 (98XJ):
    # X is correct but wrong position
    # In guess 5, X is in position 2
    if solution[2] == 'X':  # X can't be in position 2
        return False
        
    # From guess 6 (37XW):
    # Both X and W are in wrong positions
    # In guess 6, X is in position 2 and W is in position 3
    if solution[2] == 'X' or solution[3] == 'W':  # Neither can be in same positions as guess 6
        return False
        
    # From all guesses:
    # Numbers 1 or 6 must work with all number-related feedback
    # For guess 4, since 2 is in wrong position, the other number must be wrong
    if solution[0] == 1:  # Check if 1 works with all conditions
        return False  # 1 would be too small based on feedback
        
    return True

solutions = [
    [1, 2, 'X', 'W'],
    [1, 2, 'W', 'X'],
    [6, 2, 'X', 'W'],
    [6, 2, 'W', 'X']
]

for sol in solutions:
    if final_check(sol):
        print(sol)
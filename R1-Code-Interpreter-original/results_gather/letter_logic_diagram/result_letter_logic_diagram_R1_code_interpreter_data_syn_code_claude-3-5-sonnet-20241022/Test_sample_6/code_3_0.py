def is_solvable(initial_grid):
    # Check minor diagonal conflicts
    for i in range(7):
        if initial_grid[i][6-i] not in [None, ''] and initial_grid[i][6-i] != 'a':
            print(f"Conflict at row {i+1}, column {7-i}: has '{initial_grid[i][6-i]}' but needs 'a'")
            return False
    
    # Check row conflicts
    for i in range(7):
        letters = [x for x in initial_grid[i] if x not in [None, '']]
        if len(letters) != len(set(letters)):
            print(f"Duplicate letters in row {i+1}")
            return False
    
    # Check column conflicts
    for j in range(7):
        letters = [initial_grid[i][j] for i in range(7) if initial_grid[i][j] not in [None, '']]
        if len(letters) != len(set(letters)):
            print(f"Duplicate letters in column {j+1}")
            return False
    
    return True

# Initialize grid with given values
initial_grid = [
    ['c','b','a','d',None,'f',None],
    ['b',None,None,'e',None,None,None],
    [None,'d',None,'f','g',None,None],
    ['d',None,None,'g',None,None,None],
    [None,None,'g',None,None,None,None],
    [None,'g','c',None,'a',None,None],
    [None,None,'b',None,'d',None,'f']
]

print("Checking if puzzle is solvable:")
is_solvable(initial_grid)
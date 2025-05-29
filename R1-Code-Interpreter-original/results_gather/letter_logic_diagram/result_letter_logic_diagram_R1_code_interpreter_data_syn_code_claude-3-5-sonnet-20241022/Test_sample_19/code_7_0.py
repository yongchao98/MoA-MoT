def verify_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

# Initial grid
initial = [
    ['g','b','f','','c','a',''],
    ['','f','d','','a','e',''],
    ['','d','c','','','',''],
    ['','','','','','',''],
    ['','','e','','b','',''],
    ['a','e','','','','d','c'],
    ['e','','','','','c','']
]

# Create a solution starting with 'e' on diagonal
solution = [
    ['g','b','f','d','c','a','e'],  # row 0
    ['d','f','a','c','b','e','g'],  # row 1
    ['c','a','d','b','f','e','b'],  # row 2
    ['b','c','g','e','d','f','a'],  # row 3
    ['f','d','e','a','g','b','e'],  # row 4
    ['a','e','b','f','c','d','e'],  # row 5
    ['e','g','c','g','a','b','d']   # row 6
]

# Force 'e' on diagonal and maintain initial values
for i in range(7):
    # Set diagonal to 'e'
    solution[i][6-i] = 'e'

# Verify and adjust rows to maintain uniqueness
def adjust_solution():
    # First, ensure all initial values are maintained
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                solution[i][j] = initial[i][j]
    
    # Now adjust remaining positions
    letters = 'abcdefg'
    
    # Row by row adjustment
    for i in range(7):
        used = set()
        # First mark fixed positions
        for j in range(7):
            if (i, j) in [(x, 6-x) for x in range(7)] or initial[i][j] != '':
                used.add(solution[i][j])
        
        # Fill remaining positions
        available = set(letters) - used
        for j in range(7):
            if solution[i][j] not in used and (i, j) not in [(x, 6-x) for x in range(7)] and initial[i][j] == '':
                if available:
                    solution[i][j] = available.pop()
                    used.add(solution[i][j])

adjust_solution()

# Print final solution
for row in solution:
    print(','.join(row))
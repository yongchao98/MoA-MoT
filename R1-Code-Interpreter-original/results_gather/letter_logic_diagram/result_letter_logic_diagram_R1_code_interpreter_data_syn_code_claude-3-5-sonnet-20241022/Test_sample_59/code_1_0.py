def verify_prefilled(solution, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and initial[i][j] != solution[i][j]:
                return False
    return True

def verify_rows_cols(grid):
    # Check each row and column has exactly one of each letter
    letters = set('abcdefg')
    for i in range(7):
        if set(grid[i]) != letters:  # check row
            return False
        col = [grid[j][i] for j in range(7)]
        if set(col) != letters:  # check column
            return False
    return True

def create_solution():
    # Start with empty grid
    solution = [['' for _ in range(7)] for _ in range(7)]
    
    # Step 1: Fill minor diagonal with 'a' (we know this from position (0,6))
    for i in range(7):
        solution[i][6-i] = 'a'
    
    # Step 2: Fill in all pre-filled positions from initial grid
    prefilled = [
        ('0,1', 'e'), ('0,3', 'g'), ('0,4', 'c'), ('0,5', 'f'), ('0,6', 'a'),
        ('1,6', 'b'),
        ('2,3', 'f'), ('2,6', 'e'),
        ('3,3', 'a'), ('3,5', 'e'),
        ('4,3', 'b'), ('4,5', 'd'),
        ('5,0', 'f'), ('5,2', 'b'), ('5,3', 'e'), ('5,6', 'c'),
        ('6,3', 'd'), ('6,4', 'g'), ('6,6', 'f')
    ]
    
    # Fill remaining positions systematically
    solution = [
        ['b','e','d','g','c','f','a'],
        ['c','f','g','d','e','a','a'],
        ['g','d','c','f','a','b','a'],
        ['d','g','f','a','b','e','a'],
        ['e','c','a','b','f','d','a'],
        ['f','a','b','e','d','g','a'],
        ['a','b','e','d','g','c','a']
    ]
    
    return solution

# Initial grid
initial_grid = [
    ['', 'e', '', 'g', 'c', 'f', 'a'],
    ['', '', '', '', '', '', 'b'],
    ['', '', '', 'f', '', '', 'e'],
    ['', '', '', 'a', '', 'e', ''],
    ['', '', '', 'b', '', 'd', ''],
    ['f', '', 'b', 'e', '', '', 'c'],
    ['', '', '', 'd', 'g', '', 'f']
]

# Create and verify solution
solution = create_solution()

# Verify all constraints
minor_diagonal_ok = all(solution[i][6-i] == 'a' for i in range(7))
prefilled_ok = verify_prefilled(solution, initial_grid)
rows_cols_ok = verify_rows_cols(solution)

if minor_diagonal_ok and prefilled_ok and rows_cols_ok:
    for row in solution:
        print(','.join(row))
else:
    print("No valid solution found")
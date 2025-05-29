def verify_solution(grid):
    # Check rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7:  # check row
            return False
        col = [grid[j][i] for j in range(7)]
        if len(set(col)) != 7:  # check column
            return False
    # Check diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:  # all should be same
        return False
    return True

def matches_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and initial[i][j] != grid[i][j]:
                return False
    return True

# Initial grid with pre-filled positions
initial = [
    ['g','b','f','','c','a',''],
    ['','f','d','','a','e',''],
    ['','d','c','','','',''],
    ['','','','','','',''],
    ['','','e','','b','',''],
    ['a','e','','','','d','c'],
    ['e','','','','','c','']
]

# Create solution with 'e' on diagonal
solution = [
    ['g','b','f','d','c','a','e'],
    ['c','f','d','g','a','e','b'],
    ['b','d','c','f','g','e','a'],
    ['f','a','g','e','d','b','c'],
    ['d','c','e','b','b','g','e'],
    ['a','e','b','c','f','d','e'],
    ['e','g','a','d','h','c','f']
]

# Force 'e' on all diagonal positions
for i in range(7):
    solution[i][6-i] = 'e'

# Now adjust the rest to maintain row/column uniqueness
def fix_solution():
    # Known positions from initial grid
    fixed_positions = set()
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                solution[i][j] = initial[i][j]
                fixed_positions.add((i,j))
    
    # Force diagonal to be 'e'
    for i in range(7):
        solution[i][6-i] = 'e'
        fixed_positions.add((i,6-i))
    
    # Fill remaining positions ensuring uniqueness
    letters = 'abcdefg'
    
    def get_available(row, col):
        used_row = set(solution[row])
        used_col = set(solution[i][col] for i in range(7))
        return [l for l in letters if l not in used_row and l not in used_col]
    
    def solve(pos=0):
        if pos == 49:
            return True
            
        row, col = pos // 7, pos % 7
        if (row, col) in fixed_positions:
            return solve(pos + 1)
            
        available = get_available(row, col)
        for letter in available:
            solution[row][col] = letter
            if solve(pos + 1):
                return True
            solution[row][col] = ''
        return False
    
    solve()

fix_solution()

# Print solution
for row in solution:
    print(','.join(row))
def verify_solution(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

def is_valid_partial(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve(grid, fixed_positions):
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return verify_solution(grid)
        
        row, col = pos // 7, pos % 7
        
        # Skip fixed positions
        if (row, col) in fixed_positions:
            return backtrack(pos + 1)
        
        # If on diagonal, must be 'g'
        if row + col == 6:
            grid[row][col] = 'g'
            if backtrack(pos + 1):
                return True
            grid[row][col] = ''
            return False
        
        for letter in 'abcdefg':
            if is_valid_partial(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    return backtrack()

# Initial grid
initial_grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# Create working grid and track fixed positions
grid = [[''] * 7 for _ in range(7)]
fixed_positions = set()
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]
            fixed_positions.add((i, j))

# Force all diagonal positions to be 'g'
for i in range(7):
    j = 6 - i
    if (i, j) not in fixed_positions:
        grid[i][j] = 'g'
        fixed_positions.add((i, j))

if solve(grid, fixed_positions):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
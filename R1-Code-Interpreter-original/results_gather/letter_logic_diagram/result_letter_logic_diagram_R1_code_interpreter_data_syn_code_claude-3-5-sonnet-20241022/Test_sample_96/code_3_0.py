def get_diagonal_positions():
    return [(i, 6-i) for i in range(7)]

def find_fixed_diagonal_letter(grid):
    diagonal_positions = get_diagonal_positions()
    for i, j in diagonal_positions:
        if grid[i][j] != '':
            return grid[i][j]
    return 'e'  # Default if no diagonal letter is fixed

def is_valid(grid, row, col, letter, diagonal_letter):
    # If position is on diagonal, must match diagonal_letter
    if col == 6-row and letter != diagonal_letter:
        return False
        
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    return True

def solve(grid):
    # First, determine the diagonal letter and fill all diagonal positions
    diagonal_letter = find_fixed_diagonal_letter(grid)
    diagonal_positions = set((i, 6-i) for i in range(7))
    
    # Fill all diagonal positions first
    for i, j in diagonal_positions:
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip diagonal positions and pre-filled cells
        if (row, col) in diagonal_positions or grid[row][col] != '':
            return backtrack(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    return backtrack()

# Initialize grid
initial_grid = [
    ['g','','','','','',''],
    ['','','c','','','b','g'],
    ['d','','','','','g',''],
    ['c','','','','','f',''],
    ['e','','b','g','f','','c'],
    ['a','b','','f','','',''],
    ['','','f','','c','','a']
]

# Create a working copy
grid = [row[:] for row in initial_grid]

if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")
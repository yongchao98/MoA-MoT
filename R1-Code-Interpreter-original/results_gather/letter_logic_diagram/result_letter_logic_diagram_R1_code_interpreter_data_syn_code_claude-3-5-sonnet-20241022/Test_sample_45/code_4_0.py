def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_diagonal_letter(grid):
    # First check if we have any pre-filled diagonal positions
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return 'b'  # Default to 'b' if no pre-filled diagonal positions

def is_valid(grid, row, col, letter, diag_letter):
    # If this is a diagonal position, must match diagonal letter
    if row + col == 6 and letter != diag_letter:
        return False
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_grid(grid):
    diag_letter = get_diagonal_letter(grid)
    
    def backtrack():
        empty = find_empty(grid)
        if not empty:
            return True
        
        row, col = empty
        letters = [diag_letter] if row + col == 6 else 'abcdefg'
        
        for letter in letters:
            if is_valid(grid, row, col, letter, diag_letter):
                grid[row][col] = letter
                if backtrack():
                    return True
                grid[row][col] = ''
        return False
    
    return backtrack()

# Initialize grid with given values
grid = [
    ['c','','a','e','','b',''],
    ['','a','','f','','',''],
    ['','','f','b','g','',''],
    ['','','','','','',''],
    ['','','g','','','',''],
    ['b','','','d','a','e','f'],
    ['','','d','a','','f','b']
]

if solve_grid(grid):
    print_grid(grid)
else:
    print("No solution found")
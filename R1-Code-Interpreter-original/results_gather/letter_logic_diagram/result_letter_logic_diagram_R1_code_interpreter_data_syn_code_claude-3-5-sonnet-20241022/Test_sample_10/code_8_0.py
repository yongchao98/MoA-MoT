def print_solution(grid):
    for row in grid:
        print(','.join(row))

def get_prefilled_positions(grid):
    filled = []
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                filled.append((i, j, grid[i][j]))
    return filled

def is_valid_placement(grid, row, col, letter, diagonal_letter):
    # If on minor diagonal, must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
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

def solve_grid(grid, diagonal_letter):
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
        
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled positions
        if grid[row][col] != '':
            return backtrack(pos + 1)
        
        # Get available letters
        if row + col == 6:
            letters = [diagonal_letter]
        else:
            letters = [c for c in 'abcdefg']
        
        for letter in letters:
            if is_valid_placement(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    return backtrack()

# Initial grid
initial_grid = [
    ['d','','','e','b','','g'],
    ['','','e','b','','','d'],
    ['','','b','','','','c'],
    ['','','a','','d','c','f'],
    ['','a','g','d','','','e'],
    ['','','','','','',''],
    ['','','','f','','','a']
]

# Create a copy of the initial grid
grid = [row[:] for row in initial_grid]

# Try each possible letter for the diagonal
for diagonal_letter in 'abcdefg':
    # Check if this diagonal letter conflicts with existing positions
    valid = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diagonal_letter:
            valid = False
            break
    
    if not valid:
        continue
    
    # Reset grid
    grid = [row[:] for row in initial_grid]
    
    # Fill empty diagonal positions with the chosen letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    # Try to solve with this diagonal letter
    if solve_grid(grid, diagonal_letter):
        # Verify solution
        valid_solution = True
        for i in range(7):
            row_letters = [c for c in grid[i] if c != '']
            col_letters = [grid[j][i] for j in range(7) if grid[j][i] != '']
            if len(set(row_letters)) != len(row_letters) or len(set(col_letters)) != len(col_letters):
                valid_solution = False
                break
        
        if valid_solution:
            result = []
            for row in grid:
                result.append(','.join(row))
            print('\n'.join(result))
            break

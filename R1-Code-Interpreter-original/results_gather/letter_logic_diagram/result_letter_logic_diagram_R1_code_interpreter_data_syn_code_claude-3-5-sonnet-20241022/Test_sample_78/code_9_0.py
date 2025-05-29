def print_grid(grid):
    for row in grid:
        print(','.join(row))

def find_minor_diagonal_letter(grid):
    # Check what letters are already on the minor diagonal
    diagonal_letters = set()
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            diagonal_letters.add(grid[i][j])
    
    if len(diagonal_letters) == 0:
        return None
    elif len(diagonal_letters) == 1:
        return list(diagonal_letters)[0]
    else:
        return None  # Conflict in diagonal letters

def is_valid(grid, row, col, letter, diagonal_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6 and letter != diagonal_letter:
        return False
    
    return True

def solve(grid, diagonal_letter):
    def backtrack(pos=0):
        if pos == 49:
            return True
        
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells
        if grid[row][col] != '':
            return backtrack(pos + 1)
        
        # Get available letters
        used = set()
        for i in range(7):
            used.add(grid[row][i])  # Row
            used.add(grid[i][col])  # Column
        
        # If on minor diagonal, only try diagonal_letter
        if row + col == 6:
            if diagonal_letter not in used and is_valid(grid, row, col, diagonal_letter, diagonal_letter):
                grid[row][col] = diagonal_letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        else:
            # Try each available letter
            for letter in 'abcdefg':
                if letter not in used and is_valid(grid, row, col, letter, diagonal_letter):
                    grid[row][col] = letter
                    if backtrack(pos + 1):
                        return True
                    grid[row][col] = ''
        
        return False
    
    return backtrack()

# Initialize grid
initial_grid = [
    ['', 'e', '', '', '', 'g', ''],
    ['e', 'c', 'a', 'b', '', 'd', 'f'],
    ['', '', 'b', 'g', '', 'f', 'e'],
    ['', 'b', 'g', '', '', 'e', 'c'],
    ['b', '', 'd', 'f', '', 'c', ''],
    ['g', 'd', '', '', 'c', 'a', 'b'],
    ['d', '', 'e', '', 'a', 'b', '']
]

# Find what letter should be on the minor diagonal
diagonal_letter = find_minor_diagonal_letter(initial_grid)
if diagonal_letter is None:
    # Try each possible letter for the diagonal
    for letter in 'abcdefg':
        grid_copy = [row[:] for row in initial_grid]
        if solve(grid_copy, letter):
            print_grid(grid_copy)
            break
else:
    grid_copy = [row[:] for row in initial_grid]
    if solve(grid_copy, diagonal_letter):
        print_grid(grid_copy)
    else:
        print("No solution found")
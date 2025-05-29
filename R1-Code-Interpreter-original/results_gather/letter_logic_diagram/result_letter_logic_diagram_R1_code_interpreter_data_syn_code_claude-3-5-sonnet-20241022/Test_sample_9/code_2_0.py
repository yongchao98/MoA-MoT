def get_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If this is a minor diagonal position
    if row + col == 6:
        diagonal_letter = get_diagonal_letter(grid)
        if diagonal_letter is not None and letter != diagonal_letter:
            return False
    
    return True

def solve_grid(grid):
    # First determine diagonal letter if not already set
    diagonal_letter = get_diagonal_letter(grid)
    if diagonal_letter is None:
        # Try each possible letter for the entire diagonal
        for letter in 'abcdefg':
            can_use = True
            # Check if this letter can be used for all diagonal positions
            for i in range(7):
                if grid[i][6-i] != '' and grid[i][6-i] != letter:
                    can_use = False
                    break
                # Also check if letter already exists in row/column of empty diagonal positions
                if grid[i][6-i] == '':
                    for j in range(7):
                        if grid[i][j] == letter or grid[j][6-i] == letter:
                            can_use = False
                            break
            if can_use:
                # Fill all empty diagonal positions with this letter
                for i in range(7):
                    if grid[i][6-i] == '':
                        grid[i][6-i] = letter
                # Try to solve the rest
                if solve_remaining(grid):
                    return True
                # If failed, clear the diagonal
                for i in range(7):
                    if grid[i][6-i] == letter:
                        grid[i][6-i] = ''
        return False
    return solve_remaining(grid)

def solve_remaining(grid):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'c', '', '', 'b', '', ''],
    ['c', 'f', '', '', 'a', '', ''],
    ['', '', 'b', 'a', 'e', '', 'c'],
    ['', 'b', '', '', '', '', 'f'],
    ['', '', 'e', 'g', 'c', 'f', ''],
    ['', '', 'g', '', '', 'd', ''],
    ['e', '', '', '', '', 'b', '']
]

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")
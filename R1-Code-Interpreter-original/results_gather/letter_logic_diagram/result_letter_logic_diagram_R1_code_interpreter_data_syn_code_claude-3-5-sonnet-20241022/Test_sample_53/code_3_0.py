def get_minor_diagonal_letter(grid):
    letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            letters.add(grid[i][6-i])
    if len(letters) > 1:
        return None
    return next(iter(letters)) if letters else None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # If on minor diagonal, must match minor diagonal letter
    minor_letter = get_minor_diagonal_letter(grid)
    if row + col == 6 and minor_letter and letter != minor_letter:
        return False
    
    return True

def find_empty(grid):
    # First check minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then check other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    # First ensure minor diagonal has same letter
    minor_letter = get_minor_diagonal_letter(grid)
    if not minor_letter:
        # Try each letter for minor diagonal
        for letter in 'abcdefg':
            valid = True
            # Check if this letter works for all empty minor diagonal positions
            temp_grid = [row[:] for row in grid]
            for i in range(7):
                if temp_grid[i][6-i] == '':
                    if not is_valid(temp_grid, i, 6-i, letter):
                        valid = False
                        break
                    temp_grid[i][6-i] = letter
            if valid:
                # Apply changes to real grid
                for i in range(7):
                    if grid[i][6-i] == '':
                        grid[i][6-i] = letter
                if solve_remaining(grid):
                    return True
                # Backtrack
                for i in range(7):
                    if grid[i][6-i] == letter:
                        grid[i][6-i] = ''
    else:
        return solve_remaining(grid)
    return False

def solve_remaining(grid):
    empty = find_empty(grid)
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
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', 'f', ''],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")
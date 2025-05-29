def is_valid(grid, row, col, letter, diagonal_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal requirement
    if row + col == 6:  # If on diagonal
        return letter == diagonal_letter
    
    return True

def solve_with_diagonal(grid, diagonal_letter):
    def find_empty():
        # First fill diagonal positions
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i)
        # Then other positions
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None

    pos = find_empty()
    if not pos:
        return True
    
    row, col = pos
    if row + col == 6:  # Diagonal position
        if is_valid(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve_with_diagonal(grid, diagonal_letter):
                return True
            grid[row][col] = ''
        return False
    else:  # Non-diagonal position
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if solve_with_diagonal(grid, diagonal_letter):
                    return True
                grid[row][col] = ''
        return False

def check_initial_diagonal(grid):
    # Return the letter if diagonal already has a letter, None otherwise
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            if diagonal_letter is None:
                diagonal_letter = grid[i][6-i]
            elif grid[i][6-i] != diagonal_letter:
                return None  # Inconsistent diagonal
    return diagonal_letter

# Initial grid
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

print('<<<')
# Check if diagonal already has a letter
diagonal_letter = check_initial_diagonal(grid)
if diagonal_letter:
    # If diagonal letter exists, use it
    if solve_with_diagonal(grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution exists")
else:
    # Try each possible diagonal letter
    solved = False
    for d in 'abcdefg':
        grid_copy = [row[:] for row in grid]  # Make a copy of the grid
        if solve_with_diagonal(grid_copy, d):
            for row in grid_copy:
                print(','.join(row))
            solved = True
            break
    if not solved:
        print("No solution exists")
print('>>>')
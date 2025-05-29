def get_minor_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
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
    
    # If on minor diagonal, must match the diagonal letter
    if row + col == 6:
        diagonal_letter = get_minor_diagonal_letter(grid)
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve_grid(grid):
    def solve(row=0, col=0):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True

        if grid[row][col] != '':
            # For pre-filled cells on minor diagonal, verify they match
            if row + col == 6:
                diagonal_letter = get_minor_diagonal_letter(grid)
                if diagonal_letter != grid[row][col]:
                    return False
            return solve(row, col + 1)

        # Get diagonal letter if we're on the diagonal
        diagonal_letter = None
        if row + col == 6:
            diagonal_letter = get_minor_diagonal_letter(grid)

        for letter in 'abcdefg':
            # Skip if we're on diagonal and letter doesn't match
            if diagonal_letter and row + col == 6 and letter != diagonal_letter:
                continue
                
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    return solve()

# Initialize grid
initial_grid = [
    ['g', 'b', '', '', '', 'a', 'c'],
    ['b', 'e', 'f', 'd', 'a', 'c', 'g'],
    ['', '', '', '', 'c', 'g', ''],
    ['f', '', '', '', '', '', 'e'],
    ['', 'a', 'c', '', 'b', 'e', 'f'],
    ['a', '', 'g', 'b', 'e', 'f', ''],
    ['c', 'g', 'b', '', '', 'd', '']
]

# Solve the puzzle
if solve_grid(initial_grid):
    result = []
    for row in initial_grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution found")
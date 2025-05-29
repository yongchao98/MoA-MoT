def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for letter in 'abcdefg':
            if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
                continue
            grid[row][col] = letter
            if backtrack(row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    return backtrack(0, 0)

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    
    for letter in letters:
        if all(grid[i][j] != letter for i, j in minor_diagonal_positions):
            minor_diagonal_letter = letter
            break
    
    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = minor_diagonal_letter
    
    # Solve the grid
    if solve(grid, minor_diagonal_letter):
        return grid
    else:
        return None

# Initial grid
grid = [
    ['', '', 'b', '', '', 'd', 'a'],
    ['', '', 'e', '', 'd', '', 'g'],
    ['', 'e', 'c', 'd', '', '', ''],
    ['', 'c', '', '', '', '', ''],
    ['', 'd', 'a', '', 'f', '', 'e'],
    ['', '', 'g', '', '', '', ''],
    ['', '', 'f', '', 'e', 'c', 'd']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")
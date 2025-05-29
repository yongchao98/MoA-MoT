def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[row][col] = ''
                return False
    return True

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            minor_diagonal_letter = letter
            break
    
    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
    
    # Solve the grid using backtracking
    solve(grid, minor_diagonal_letter)
    return grid

# Initial grid setup
grid = [
    ['c', 'b', '', 'a', 'f', 'g', ''],
    ['b', 'e', 'a', 'f', '', '', ''],
    ['', '', '', 'g', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['', 'g', '', 'c', 'b', '', ''],
    ['g', 'd', '', '', '', '', 'f'],
    ['', '', 'b', '', 'a', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    # Try to fill the grid using backtracking
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
    possible_letters = set('abcdefg')
    for i, j in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        if grid[i][j] != '':
            possible_letters.intersection_update(grid[i][j])
    
    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal
    for i, j in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        grid[i][j] = minor_diagonal_letter
    
    # Solve the grid
    solve(grid, minor_diagonal_letter)
    
    return grid

# Initial grid
grid = [
    ['b', '', 'e', '', '', '', ''],
    ['', '', '', '', 'f', '', ''],
    ['', '', '', 'f', 'g', '', ''],
    ['c', '', 'f', 'g', 'b', 'd', ''],
    ['', 'f', '', 'b', 'd', '', 'c'],
    ['', 'g', '', 'd', '', '', 'a'],
    ['g', '', '', '', '', '', 'f']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
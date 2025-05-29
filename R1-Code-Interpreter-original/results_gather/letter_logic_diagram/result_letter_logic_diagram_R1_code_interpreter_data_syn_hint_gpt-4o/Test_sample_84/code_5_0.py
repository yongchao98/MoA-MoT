def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, diagonal_letter, row, col + 1)
    
    # If on the minor diagonal, place the diagonal letter
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Try placing each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_diagonal_letters = set('abcdefg')
    
    for r, c in diagonal_positions:
        if grid[r][c] != '':
            possible_diagonal_letters.intersection_update(grid[r][c])
    
    diagonal_letter = possible_diagonal_letters.pop()
    
    # Solve the grid
    solve(grid, diagonal_letter)
    return grid

# Initial grid
grid = [
    ['', '', 'f', '', '', 'e', 'c'],
    ['d', 'f', 'b', 'g', 'e', 'c', ''],
    ['f', 'b', 'g', 'e', '', 'a', ''],
    ['', 'g', '', '', 'a', 'd', ''],
    ['', 'e', '', '', 'd', '', 'b'],
    ['e', '', 'a', 'd', 'f', 'b', 'g'],
    ['', '', 'd', 'f', '', 'g', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
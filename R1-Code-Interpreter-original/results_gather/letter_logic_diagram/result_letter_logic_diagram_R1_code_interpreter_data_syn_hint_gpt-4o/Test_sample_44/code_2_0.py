def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_diagonal_letters = set('abcdefg')
    
    for r, c in minor_diagonal_positions:
        if grid[r][c] != '':
            possible_diagonal_letters.intersection_update(grid[r][c])
    
    diagonal_letter = possible_diagonal_letters.pop()
    
    for r, c in minor_diagonal_positions:
        grid[r][c] = diagonal_letter
    
    if solve(grid):
        return grid
    else:
        return None

# Initial grid
grid = [
    ['b', '', '', '', '', 'a', 'f'],
    ['', '', '', '', '', '', 'b'],
    ['', 'e', 'd', 'a', '', '', ''],
    ['e', '', '', 'f', 'b', 'g', ''],
    ['', 'a', '', 'b', 'g', '', 'e'],
    ['', '', 'b', 'g', '', '', 'd'],
    ['', '', '', '', '', 'd', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid in the required format
if filled_grid:
    print("<<<")
    for row in filled_grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution found")
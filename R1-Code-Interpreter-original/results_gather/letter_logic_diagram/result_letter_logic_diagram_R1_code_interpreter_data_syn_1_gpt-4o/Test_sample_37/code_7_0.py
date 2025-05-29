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

def solve_puzzle(grid):
    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            diagonal_letter = letter
            break
    
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter
    
    # Solve the rest of the grid
    solve(grid)
    return grid

# Initial grid
grid = [
    ['c', '', 'g', '', 'f', 'e', ''],
    ['a', 'g', '', '', '', 'd', ''],
    ['', '', '', '', 'd', '', ''],
    ['', '', '', 'd', '', 'a', 'g'],
    ['', 'e', '', 'c', 'a', 'g', 'b'],
    ['e', '', '', 'a', 'g', 'b', 'f'],
    ['d', '', 'a', '', 'b', '', '']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))
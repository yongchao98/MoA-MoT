def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the current row
    if letter in grid[row]:
        return False
    # Check if the letter can be placed in the current column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid):
    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        # Place the letter on the minor diagonal
        for i in range(7):
            grid[i][6-i] = letter
        
        # Try to fill the rest of the grid
        if fill_grid(grid, 0, 0):
            return grid
    return None

def fill_grid(grid, row, col):
    if row == 7:
        return True
    if col == 7:
        return fill_grid(grid, row + 1, 0)
    if grid[row][col] != '':
        return fill_grid(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if fill_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid with pre-filled values
grid = [
    ['d', 'f', 'e', 'b', 'a', 'g', ''],
    ['f', '', 'b', 'a', 'g', '', ''],
    ['e', 'b', '', '', 'c', 'd', ''],
    ['b', '', 'g', 'c', '', 'f', ''],
    ['b', 'a', '', 'c', 'd', '', ''],
    ['a', '', '', 'd', 'f', 'e', 'b'],
    ['', '', '', '', 'e', '', 'a']
]

# Solve the grid
solution = solve(grid)

# Check if a solution was found and print it
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")
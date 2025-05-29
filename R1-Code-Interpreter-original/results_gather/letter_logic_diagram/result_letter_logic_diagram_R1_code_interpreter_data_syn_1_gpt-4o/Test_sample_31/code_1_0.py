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
    
    # Determine the letter for the minor diagonal
    if row + col == 6:
        for letter in 'abcdefg':
            if all(grid[i][6-i] == letter or grid[i][6-i] == '' for i in range(7)):
                if is_valid(grid, row, col, letter):
                    grid[row][col] = letter
                    if solve(grid, row, col + 1):
                        return True
                    grid[row][col] = ''
        return False
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['e', 'a', '', '', '', 'c', 'f'],
    ['', 'b', '', '', 'c', 'f', ''],
    ['b', 'g', 'd', 'c', 'f', '', 'a'],
    ['g', 'd', '', 'f', 'e', '', 'b'],
    ['d', '', '', 'e', '', '', 'g'],
    ['', 'f', 'e', '', 'b', '', 'd'],
    ['', '', 'a', '', 'g', 'd', 'c']
]

# Solve the puzzle
solve(grid)

# Print the solved grid
for row in grid:
    print(','.join(row))
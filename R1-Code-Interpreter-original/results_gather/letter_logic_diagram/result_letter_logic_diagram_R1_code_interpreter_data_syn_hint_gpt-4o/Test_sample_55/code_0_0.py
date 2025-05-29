def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    # Check if the letter is valid for the minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] not in ('', letter):
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

# Initial grid setup
grid = [
    ['a', 'g', '', '', 'f', '', 'b'],
    ['', 'd', '', '', 'e', '', 'a'],
    ['d', '', 'f', 'e', '', 'a', 'g'],
    ['c', 'f', 'e', 'b', 'a', 'g', ''],
    ['', '', 'b', '', 'g', 'd', 'c'],
    ['e', 'b', 'a', 'g', 'd', 'c', ''],
    ['b', '', 'g', 'd', 'c', 'f', 'e']
]

# Solve the puzzle
solve(grid)

# Print the solution
for row in grid:
    print(','.join(row))
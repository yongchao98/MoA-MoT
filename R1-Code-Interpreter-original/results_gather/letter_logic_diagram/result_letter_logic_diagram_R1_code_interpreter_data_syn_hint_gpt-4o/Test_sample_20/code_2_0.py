def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)
    
    # If the cell is on the minor diagonal, it must be the minor_diagonal_letter
    if col == 6 - row:
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Try placing each letter in the current cell
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid setup
grid = [
    ['b', '', '', 'g', '', 'c', 'e'],
    ['', 'd', 'g', 'a', '', 'e', ''],
    ['d', '', 'a', 'c', '', 'b', 'f'],
    ['g', '', '', 'e', '', '', 'd'],
    ['', '', 'e', 'b', '', '', ''],
    ['', 'e', '', 'f', 'd', 'g', 'a'],
    ['e', 'b', 'f', 'd', 'g', '', 'c']
]

# Choose a letter for the minor diagonal
minor_diagonal_letter = 'a'

# Solve the puzzle
solve(grid, minor_diagonal_letter)

# Print the filled grid
for row in grid:
    print(','.join(row))
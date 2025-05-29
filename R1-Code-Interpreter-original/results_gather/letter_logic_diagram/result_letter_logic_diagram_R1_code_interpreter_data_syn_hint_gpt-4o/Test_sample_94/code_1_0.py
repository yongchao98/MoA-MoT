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
    
    # If the current position is on the minor diagonal, place the minor diagonal letter
    if col == 6 - row:
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Try placing each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid setup
grid = [
    ['a', 'f', 'g', 'c', 'e', 'd', 'b'],
    ['f', '', 'c', '', '', 'b', 'a'],
    ['g', 'c', 'e', '', '', 'a', 'f'],
    ['c', '', '', '', '', 'f', ''],
    ['e', 'd', '', 'a', 'f', '', ''],
    ['d', 'b', 'a', '', 'g', '', 'e'],
    ['b', '', '', '', '', '', 'd']
]

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    if solve(grid, letter):
        break

# Print the filled grid
for row in grid:
    print(','.join(row))
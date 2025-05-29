def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    # Check if the letter is valid for the minor diagonal
    if row + col == 6 and letter != diagonal_letter:
        return False
    return True

def solve(grid, row=0, col=0, diagonal_letter=None):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0, diagonal_letter)
    if grid[row][col] != '':
        return solve(grid, row, col + 1, diagonal_letter)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1, diagonal_letter):
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

# Determine the diagonal letter
# We can choose 'a' as it appears frequently and is already on the diagonal
diagonal_letter = 'a'

# Solve the puzzle
solve(grid, diagonal_letter=diagonal_letter)

# Print the solution
for row in grid:
    print(','.join(row))
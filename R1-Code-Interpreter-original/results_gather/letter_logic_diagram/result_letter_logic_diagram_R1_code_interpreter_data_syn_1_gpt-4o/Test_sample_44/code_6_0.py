def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for letter in 'abcdefg':
            if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
                continue
            grid[row][col] = letter
            if backtrack(row, col + 1):
                return True
            grid[row][col] = ''
        return False

    backtrack(0, 0)
    return grid

# Initial grid setup
grid = [
    ['b', '', '', '', '', 'a', 'f'],
    ['', '', '', '', '', '', 'b'],
    ['', 'e', 'd', 'a', '', '', ''],
    ['e', '', '', 'f', 'b', 'g', ''],
    ['', 'a', '', 'b', 'g', '', 'e'],
    ['', '', 'b', 'g', '', '', 'd'],
    ['', '', '', '', '', 'd', '']
]

# Determine the minor diagonal letter
# We need to choose a letter that is not already on the minor diagonal
# Check the current minor diagonal positions: (0,6), (1,5), (2,4), (3,3), (4,2), (5,1), (6,0)
# Current letters on the minor diagonal: f, '', '', f, '', '', ''
# Choose a letter that is not 'f' and not already present in these positions
minor_diagonal_letter = 'c'  # Example choice, adjust if needed
solved_grid = solve(grid, minor_diagonal_letter)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))
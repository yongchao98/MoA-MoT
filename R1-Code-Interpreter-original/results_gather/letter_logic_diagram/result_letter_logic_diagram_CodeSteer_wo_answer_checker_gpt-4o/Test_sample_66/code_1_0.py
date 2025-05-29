def is_valid(grid, row, col, letter, rows, cols, diag_letter):
    # Check if the letter is already in the row or column
    if letter in rows[row] or letter in cols[col]:
        return False
    # Check if the letter is valid for the minor diagonal
    if row + col == 6 and letter != diag_letter:
        return False
    return True

def solve(grid, rows, cols, diag_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, rows, cols, diag_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, rows, cols, diag_letter, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, rows, cols, diag_letter):
            grid[row][col] = letter
            rows[row].add(letter)
            cols[col].add(letter)
            if solve(grid, rows, cols, diag_letter, row, col + 1):
                return True
            grid[row][col] = ''
            rows[row].remove(letter)
            cols[col].remove(letter)
    return False

# Initial grid setup
grid = [
    ['', 'c', '', 'e', '', '', ''],
    ['c', '', 'e', '', '', '', 'd'],
    ['a', '', 'b', 'f', '', '', 'c'],
    ['', '', '', '', '', 'c', 'a'],
    ['b', '', 'g', 'd', '', '', ''],
    ['', 'g', '', '', '', 'e', 'b'],
    ['', '', 'c', '', 'e', '', '']
]

# Determine the letter for the minor diagonal
diag_letter = None
for i in range(7):
    if grid[i][6-i] != '':
        diag_letter = grid[i][6-i]
        break
if diag_letter is None:
    diag_letter = 'a'  # Choose 'a' if no pre-filled diagonal letter

# Initialize sets for rows and columns
rows = [set(row) - {''} for row in grid]
cols = [set(grid[i][j] for i in range(7)) - {''} for j in range(7)]

# Solve the puzzle
solve(grid, rows, cols, diag_letter)

# Print the solution
for row in grid:
    print(','.join(row))
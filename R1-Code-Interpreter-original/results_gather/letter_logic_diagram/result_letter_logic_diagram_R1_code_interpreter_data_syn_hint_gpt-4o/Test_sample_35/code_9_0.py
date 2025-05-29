def is_valid(grid, row, col, letter, minor_diagonal_letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    # Check the minor diagonal constraint
    if row + col == 6 and letter != minor_diagonal_letter:
        return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter, minor_diagonal_letter):
                        grid[row][col] = letter
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[row][col] = ''
                return False
    return True

# Initial grid setup
grid = [
    ['', 'a', '', 'g', 'b', '', ''],
    ['a', 'd', 'g', '', '', 'c', ''],
    ['d', '', '', 'e', 'c', '', 'a'],
    ['', 'b', 'e', '', 'f', 'a', 'd'],
    ['b', 'e', 'c', '', '', 'd', 'g'],
    ['', '', '', '', 'd', 'g', ''],
    ['', '', 'a', '', '', '', '']
]

# Try each letter for the minor diagonal
for minor_diagonal_letter in 'abcdefg':
    # Check if the minor diagonal letter is valid with pre-filled cells
    valid = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != minor_diagonal_letter:
            valid = False
            break
    if not valid:
        continue

    # Set the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    if solve(grid, minor_diagonal_letter):
        break

# Print the solved grid
for row in grid:
    print(','.join(row))
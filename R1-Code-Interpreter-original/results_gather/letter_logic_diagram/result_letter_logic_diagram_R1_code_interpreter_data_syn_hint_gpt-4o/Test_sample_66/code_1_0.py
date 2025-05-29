def is_valid(grid, row, col, letter, row_sets, col_sets):
    # Check if the letter can be placed at grid[row][col]
    return letter not in row_sets[row] and letter not in col_sets[col]

def solve(grid, row_sets, col_sets, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row_sets, col_sets, minor_diagonal_letter, row + 1, 0)

    if grid[row][col] != '':
        return solve(grid, row_sets, col_sets, minor_diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row, col) in minor_diagonal_indices and letter != minor_diagonal_letter:
            continue

        if is_valid(grid, row, col, letter, row_sets, col_sets):
            grid[row][col] = letter
            row_sets[row].add(letter)
            col_sets[col].add(letter)

            if solve(grid, row_sets, col_sets, minor_diagonal_letter, row, col + 1):
                return True

            grid[row][col] = ''
            row_sets[row].remove(letter)
            col_sets[col].remove(letter)

    return False

# Initial grid
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
minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
possible_letters = set('abcdefg')

for r, c in minor_diagonal_indices:
    if grid[r][c] != '':
        possible_letters.intersection_update(grid[r][c])

minor_diagonal_letter = possible_letters.pop()

# Initialize row and column sets
row_sets = [set(row) - {''} for row in grid]
col_sets = [set(grid[r][c] for r in range(7)) - {''} for c in range(7)]

# Solve the grid
solve(grid, row_sets, col_sets, minor_diagonal_letter)

# Print the filled grid
for row in grid:
    print(','.join(row))
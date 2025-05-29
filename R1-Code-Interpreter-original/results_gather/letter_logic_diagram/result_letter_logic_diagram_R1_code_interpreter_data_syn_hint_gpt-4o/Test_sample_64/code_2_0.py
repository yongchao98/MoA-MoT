def is_valid(grid, row, col, letter, row_sets, col_sets):
    # Check if the letter can be placed at grid[row][col]
    return letter in row_sets[row] and letter in col_sets[col]

def solve(grid, minor_diagonal_letter, row_sets, col_sets, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_letter, row_sets, col_sets, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row_sets, col_sets, row, col + 1)
    
    for letter in 'abcdefg':
        if (row, col) in minor_diagonal_indices and letter != minor_diagonal_letter:
            continue
        if is_valid(grid, row, col, letter, row_sets, col_sets):
            grid[row][col] = letter
            row_sets[row].remove(letter)
            col_sets[col].remove(letter)
            if solve(grid, minor_diagonal_letter, row_sets, col_sets, row, col + 1):
                return True
            grid[row][col] = ''
            row_sets[row].add(letter)
            col_sets[col].add(letter)
    return False

# Initial grid with pre-filled values
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Determine the letter for the minor diagonal
minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
minor_diagonal_letter = None

for i, j in minor_diagonal_indices:
    if initial_grid[i][j] != '':
        minor_diagonal_letter = initial_grid[i][j]
        break

if minor_diagonal_letter is None:
    for letter in 'abcdefg':
        if all(initial_grid[i][j] != letter for i, j in minor_diagonal_indices):
            minor_diagonal_letter = letter
            break

# Initialize row and column sets
row_sets = [set('abcdefg') - set(row) for row in initial_grid]
col_sets = [set('abcdefg') - {initial_grid[i][j] for i in range(7)} for j in range(7)]

# Solve the puzzle
solve(initial_grid, minor_diagonal_letter, row_sets, col_sets)

# Format the output
solution = '\n'.join(','.join(row) for row in initial_grid)
print(solution)
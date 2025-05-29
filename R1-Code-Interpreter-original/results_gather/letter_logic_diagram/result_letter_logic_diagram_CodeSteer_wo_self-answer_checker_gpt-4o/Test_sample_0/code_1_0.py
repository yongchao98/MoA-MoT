def is_valid(grid, row, col, letter, row_sets, col_sets):
    # Check if the letter can be placed at grid[row][col]
    if letter in row_sets[row] or letter in col_sets[col]:
        return False
    return True

def find_unassigned_location(grid, row_sets, col_sets):
    # Find the unassigned location with the fewest possible options
    min_options = float('inf')
    best_row, best_col = -1, -1
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                options = 0
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter, row_sets, col_sets):
                        options += 1
                if options < min_options:
                    min_options = options
                    best_row, best_col = row, col
    return best_row, best_col

def solve(grid, diagonal_letter, row_sets, col_sets):
    row, col = find_unassigned_location(grid, row_sets, col_sets)
    if row == -1 and col == -1:
        return True

    for letter in 'abcdefg':
        if (row + col == 6 and letter != diagonal_letter) or not is_valid(grid, row, col, letter, row_sets, col_sets):
            continue
        grid[row][col] = letter
        row_sets[row].add(letter)
        col_sets[col].add(letter)
        print(f"Placing {letter} at ({row}, {col})")  # Debugging output
        if solve(grid, diagonal_letter, row_sets, col_sets):
            return True
        print(f"Backtracking from ({row}, {col})")  # Debugging output
        grid[row][col] = ''
        row_sets[row].remove(letter)
        col_sets[col].remove(letter)

    return False

def solve_puzzle():
    grid = [
        ['', '', 'e', '', '', '', ''],
        ['g', 'e', 'f', '', '', '', 'c'],
        ['', '', '', 'a', '', '', ''],
        ['', '', '', '', '', '', 'e'],
        ['', 'a', 'd', 'c', '', '', ''],
        ['', 'd', 'c', 'g', '', '', 'b'],
        ['d', 'c', 'g', '', 'f', 'b', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_indices if grid[i][j]}
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Initialize sets to track available letters for each row and column
    row_sets = [set(row) for row in grid]
    col_sets = [set(grid[i][j] for i in range(7)) for j in range(7)]

    # Solve the puzzle using backtracking
    if solve(grid, diagonal_letter, row_sets, col_sets):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()
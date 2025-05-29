def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check if the letter is valid for the diagonal
    if (row + col == 6) and (diagonal_letter is not None) and (letter != diagonal_letter):
        return False
    return True

def solve(grid, row, col, diagonal_letter):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0, diagonal_letter)  # Move to the next row

    if grid[row][col] is not None:
        return solve(grid, row, col + 1, diagonal_letter)  # Skip pre-filled cells

    # Determine possible letters for this cell
    possible_letters = set('abcdefg') - set(grid[row]) - set(grid[i][col] for i in range(7))

    for letter in possible_letters:
        # Dynamically choose the diagonal letter if needed
        if row + col == 6 and diagonal_letter is None:
            diagonal_letter = letter

        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            print(f"Placing {letter} at ({row}, {col}), Diagonal Letter: {diagonal_letter}")
            if solve(grid, row, col + 1, diagonal_letter):
                return True
            grid[row][col] = None  # Backtrack
            print(f"Backtracking from ({row}, {col})")

    return False

def solve_puzzle():
    grid = [
        [None, 'g', 'e', None, 'b', None, None],
        ['g', 'e', 'f', 'b', None, None, None],
        [None, None, None, None, None, None, 'g'],
        ['f', None, 'd', None, None, 'g', 'e'],
        [None, 'd', None, None, None, None, None],
        ['d', 'c', None, None, None, None, 'b'],
        ['c', None, 'g', None, None, 'b', None]
    ]

    if solve(grid, 0, 0, None):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()
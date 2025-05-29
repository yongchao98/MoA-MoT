def solve_puzzle(grid):
    def can_place(letter, row, col):
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    def forward_checking():
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    possible_letters = rows_missing[r] & cols_missing[c]
                    if not possible_letters:
                        return False
        return True

    def fill_grid():
        # Find the first empty cell
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    possible_letters = rows_missing[r] & cols_missing[c]
                    for letter in possible_letters:
                        grid[r][c] = letter
                        rows_missing[r].remove(letter)
                        cols_missing[c].remove(letter)

                        if forward_checking() and fill_grid():
                            return True

                        # Backtrack
                        grid[r][c] = ''
                        rows_missing[r].add(letter)
                        cols_missing[c].add(letter)

                    return False
        return True

    # Find the missing letters for each row and column
    rows_missing = [set('abcdefg') - set(row) for row in grid]
    cols_missing = [set('abcdefg') - set(grid[r][c] for r in range(7)) for c in range(7)]

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if all(can_place(letter, r, 6-r) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
                rows_missing[r].discard(letter)
                cols_missing[6-r].discard(letter)

            # Try to fill the rest of the grid
            if fill_grid():
                return grid

            # Reset the grid if unsuccessful
            for r in range(7):
                grid[r][6-r] = ''
                rows_missing[r].add(letter)
                cols_missing[6-r].add(letter)

    return None

# Initial grid
grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    print("<<<")
    for row in solution:
        print(','.join(row))
    print(">>>")
else:
    print("No solution found.")
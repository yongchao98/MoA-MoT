def solve_puzzle(grid):
    # Helper function to find the missing letter in a row or column
    def find_missing(letters):
        all_letters = set('abcdefg')
        return list(all_letters - set(letters))[0]

    # Fill the minor diagonal with the same letter
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            minor_diagonal_letter = grid[i][6-i]
            break

    if minor_diagonal_letter is None:
        # If no letter is pre-filled on the minor diagonal, choose one that can fit
        for letter in 'abcdefg':
            if all(letter not in grid[i] and letter not in [grid[j][6-i] for j in range(7)] for i in range(7)):
                minor_diagonal_letter = letter
                break

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find the missing letter for this row
                row_letters = [grid[i][k] for k in range(7) if grid[i][k] != '']
                col_letters = [grid[k][j] for k in range(7) if grid[k][j] != '']
                missing_letter = find_missing(row_letters + col_letters)
                grid[i][j] = missing_letter

    return grid

# Initial grid setup
grid = [
    ['f', '', '', '', 'e', 'a', 'g'],
    ['c', 'd', 'b', '', 'a', 'g', 'f'],
    ['', 'b', '', 'a', 'g', 'f', 'c'],
    ['b', 'e', '', 'g', 'f', 'c', 'd'],
    ['', 'a', 'g', 'f', 'c', 'd', 'b'],
    ['', 'g', 'f', '', '', 'b', 'e'],
    ['g', '', '', '', '', 'e', 'a']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))
def solve_puzzle(grid):
    from itertools import permutations

    # Helper function to check if the grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Extract the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    minor_diagonal = [grid[i][j] for i, j in minor_diagonal_indices if grid[i][j] != '']

    # If the minor diagonal is not consistent, return None
    if len(set(minor_diagonal)) > 1:
        return None

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = minor_diagonal[0] if minor_diagonal else None

    # Fill the grid
    letters = set('abcdefg')
    for i in range(7):
        row_letters = set(grid[i])
        missing_letters = letters - row_letters
        for perm in permutations(missing_letters):
            new_row = list(grid[i])
            perm_index = 0
            for j in range(7):
                if new_row[j] == '':
                    new_row[j] = perm[perm_index]
                    perm_index += 1
            grid[i] = new_row

    # Check if the grid is valid
    if is_valid(grid):
        # Ensure the minor diagonal is consistent
        for i, j in minor_diagonal_indices:
            grid[i][j] = minor_diagonal_letter

        # Format the output
        result = '\n'.join(','.join(row) for row in grid)
        return result
    else:
        return None

# Initial grid
grid = [
    ['c', 'b', '', 'a', 'f', 'g', ''],
    ['b', 'e', 'a', 'f', '', '', ''],
    ['', '', '', 'g', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['', 'g', '', 'c', 'b', '', ''],
    ['g', 'd', '', '', '', '', 'f'],
    ['', '', 'b', '', 'a', '', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)
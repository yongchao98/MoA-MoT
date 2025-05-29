def solve_puzzle(grid):
    import copy

    def is_valid(grid, row, col, char):
        # Check if char is not in the current row and column
        for i in range(7):
            if grid[row][i] == char or grid[i][col] == char:
                return False
        return True

    def solve(grid, minor_diagonal_letter):
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for char in 'abcdefg':
                        if is_valid(grid, row, col, char):
                            grid[row][col] = char
                            if (row, col) in minor_diagonal_indices:
                                if char != minor_diagonal_letter:
                                    grid[row][col] = ''
                                    continue
                            if solve(grid, minor_diagonal_letter):
                                return True
                            grid[row][col] = ''
                    return False
        return True

    # Extract the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    minor_diagonal = [grid[i][j] for i, j in minor_diagonal_indices if grid[i][j] != '']

    # Determine the letter for the minor diagonal
    if len(set(minor_diagonal)) > 1:
        return None
    minor_diagonal_letter = minor_diagonal[0] if minor_diagonal else 'a'  # Default to 'a' if empty

    # Make a deep copy of the grid to avoid modifying the original
    grid_copy = copy.deepcopy(grid)

    # Solve the puzzle
    if solve(grid_copy, minor_diagonal_letter):
        return '\n'.join(','.join(row) for row in grid_copy)
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
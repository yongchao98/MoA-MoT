def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, minor_diagonal_letter):
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in 'abcdefg':
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if (row, col) in minor_diagonal_indices:
                                if letter != minor_diagonal_letter:
                                    grid[row][col] = ''
                                    continue
                            if solve(grid, minor_diagonal_letter):
                                return True
                            grid[row][col] = ''
                    return False
        return True

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    minor_diagonal_letter = None
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            minor_diagonal_letter = grid[i][j]
            break

    if not minor_diagonal_letter:
        # If no pre-filled minor diagonal letter, choose one that fits
        for letter in 'abcdefg':
            if all(letter not in grid[i] and letter not in [grid[j][i] for j in range(7)] for i in range(7)):
                minor_diagonal_letter = letter
                break

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter

    # Solve the puzzle
    solve(grid, minor_diagonal_letter)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid setup
grid = [
    ['', 'd', '', '', 'b', 'e', 'g'],
    ['', '', '', 'b', '', 'g', ''],
    ['f', '', '', 'e', 'g', '', ''],
    ['a', '', '', '', '', 'd', ''],
    ['', '', 'g', '', '', '', 'a'],
    ['', 'g', '', 'd', '', '', ''],
    ['', 'c', 'd', '', 'a', 'b', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)
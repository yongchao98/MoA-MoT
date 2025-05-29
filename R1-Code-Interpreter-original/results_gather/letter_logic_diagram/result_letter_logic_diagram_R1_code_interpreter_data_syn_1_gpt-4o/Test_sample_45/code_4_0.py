def solve_puzzle(grid):
    # Helper function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        return True

    # Find the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = set('abcdefg')
    for r, c in minor_diagonal_positions:
        if grid[r][c] != '':
            possible_letters.intersection_update(grid[r][c])
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if can_place(letter, row, col):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    return grid

# Initial grid setup
grid = [
    ['c', '', 'a', 'e', '', 'b', ''],
    ['', 'a', '', 'f', '', '', ''],
    ['', '', 'f', 'b', 'g', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'g', '', '', '', ''],
    ['b', '', '', 'd', 'a', 'e', 'f'],
    ['', '', 'd', 'a', '', 'f', 'b']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))
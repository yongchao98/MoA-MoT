def solve_puzzle(grid):
    def is_valid(row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        if letter in grid[row]:  # Check row
            return False
        if any(grid[r][col] == letter for r in range(7)):  # Check column
            return False
        return True

    def solve():
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in 'abcdefg':
                        if is_valid(row, col, letter):
                            grid[row][col] = letter
                            if solve():
                                return True
                            grid[row][col] = ''  # Backtrack
                    return False
        return True

    # Determine the diagonal letter
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for r, c in diagonal_indices:
        if grid[r][c] != '':
            diagonal_letter = grid[r][c]
            break

    if not diagonal_letter:
        for letter in 'abcdefg':
            if all(is_valid(r, c, letter) for r, c in diagonal_indices):
                diagonal_letter = letter
                break

    # Fill the diagonal with the chosen letter
    for r, c in diagonal_indices:
        grid[r][c] = diagonal_letter

    # Solve the puzzle
    solve()

    # Print the solution
    for row in grid:
        print(','.join(row))

# Initial grid with empty cells as ''
initial_grid = [
    ['a', 'b', '', 'g', 'e', 'c', 'd'],
    ['', 'f', 'g', 'e', '', 'd', ''],
    ['f', 'g', 'e', '', '', 'a', ''],
    ['', 'e', 'c', 'd', '', '', 'f'],
    ['e', '', 'd', '', '', 'f', ''],
    ['c', '', 'a', 'b', '', 'g', ''],
    ['', '', 'b', '', 'g', 'e', 'c']
]

solve_puzzle(initial_grid)
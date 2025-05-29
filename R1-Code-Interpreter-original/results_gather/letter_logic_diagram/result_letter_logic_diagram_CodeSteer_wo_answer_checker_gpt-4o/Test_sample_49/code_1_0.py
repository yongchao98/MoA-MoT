def solve_puzzle(grid):
    import copy

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
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            minor_diagonal_letter = grid[i][6-i]
            break

    # If no pre-filled letter on the minor diagonal, choose one
    if not minor_diagonal_letter:
        # Choose a letter that is not already over-represented
        letter_counts = {chr(c): 0 for c in range(ord('a'), ord('h'))}
        for row in grid:
            for letter in row:
                if letter:
                    letter_counts[letter] += 1
        minor_diagonal_letter = min(letter_counts, key=letter_counts.get)

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(row):
        if row == 7:
            return True
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if can_place(letter, row, col):
                        grid[row][col] = letter
                        if backtrack(row):
                            return True
                        grid[row][col] = ''
                return False
        return backtrack(row + 1)

    # Start backtracking from the first row
    backtrack(0)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid with empty cells as ''
initial_grid = [
    ['', '', '', '', 'd', '', 'e'],
    ['c', 'f', '', '', '', 'g', ''],
    ['f', 'b', '', '', '', 'a', 'c'],
    ['', 'd', '', 'e', '', 'c', ''],
    ['', 'g', '', 'a', '', '', 'b'],
    ['', '', 'a', 'c', '', '', ''],
    ['', '', '', '', 'b', 'd', '']
]

# Solve the puzzle
solution = solve_puzzle(initial_grid)
print(solution)
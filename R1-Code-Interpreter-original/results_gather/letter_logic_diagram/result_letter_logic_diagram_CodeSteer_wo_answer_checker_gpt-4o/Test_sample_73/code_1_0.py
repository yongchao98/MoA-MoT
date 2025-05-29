def solve_puzzle(grid):
    # Helper function to check if a grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Find the letter for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][j] in ('', letter) for i, j in diagonal_indices):
            diagonal_letter = letter
            break

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Initialize sets for available letters in each row and column
    rows = [set('abcdefg') - set(row) for row in grid]
    cols = [set('abcdefg') - set(grid[i][j] for i in range(7)) for j in range(7)]

    # Backtracking function
    def backtrack(r, c):
        if r == 7:
            return is_valid(grid)
        if c == 7:
            return backtrack(r + 1, 0)
        if grid[r][c] != '':
            return backtrack(r, c + 1)

        for letter in rows[r] & cols[c]:
            grid[r][c] = letter
            rows[r].remove(letter)
            cols[c].remove(letter)

            if backtrack(r, c + 1):
                return True

            # Undo the choice
            grid[r][c] = ''
            rows[r].add(letter)
            cols[c].add(letter)

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid
grid = [
    ['f', 'e', 'g', 'a', '', '', ''],
    ['e', 'g', 'a', '', '', '', 'f'],
    ['g', '', '', '', 'b', 'f', 'e'],
    ['', 'd', '', '', 'f', 'e', 'g'],
    ['d', '', 'b', '', 'e', '', 'a'],
    ['', 'b', 'f', 'e', '', 'a', 'd'],
    ['', '', 'e', '', 'a', 'd', 'c']
]

# Solve the puzzle
print(solve_puzzle(grid))
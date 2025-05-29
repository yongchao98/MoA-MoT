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

    # Find the letter for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = [grid[i][j] for i, j in diagonal_indices if grid[i][j] != '']
    diagonal_letter = next(letter for letter in 'abcdefg' if letter not in diagonal_letters)

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Available letters for each row and column
    row_available = [set('abcdefg') - set(row) for row in grid]
    col_available = [set('abcdefg') - set(grid[row][col] for row in range(7)) for col in range(7)]

    # Try to fill the rest of the grid
    def backtrack(row):
        if row == 7:
            return is_valid(grid)
        for perm in permutations('abcdefg'):
            if all(grid[row][col] in ('', perm[col]) for col in range(7)):
                original_row = grid[row][:]
                grid[row] = [perm[col] if grid[row][col] == '' else grid[row][col] for col in range(7)]
                if backtrack(row + 1):
                    return True
                grid[row] = original_row
        return False

    # Start backtracking from the first row
    backtrack(0)

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid with empty cells as ''
grid = [
    ['c', 'b', 'a', 'd', '', 'f', ''],
    ['b', '', '', 'e', '', '', ''],
    ['', 'd', '', 'f', 'g', '', ''],
    ['d', '', '', 'g', '', '', ''],
    ['', '', 'g', '', '', '', ''],
    ['', 'g', 'c', '', 'a', '', ''],
    ['', '', 'b', '', 'd', '', 'f']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print("<<<")
print(solution)
print(">>>")
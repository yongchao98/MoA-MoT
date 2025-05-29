def solve_puzzle(grid):
    from collections import Counter

    # Convert the grid into a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(i, 6 - i) for i in range(7)]
    minor_diagonal_letters = [grid[i][6 - i] for i in range(7) if grid[i][6 - i]]
    letter_counts = Counter(minor_diagonal_letters)
    all_letters = set('abcdefg')
    diagonal_letter = min(all_letters - set(minor_diagonal_letters), key=lambda x: letter_counts.get(x, 0))

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if (row, col) in minor_diagonal_indices or grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in all_letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Start solving the puzzle
    solve(grid)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Given puzzle
puzzle = """
c,b,a,d,,f,
b,,,e,,,
,d,,f,g,,
d,,,g,,,
,,g,,,,
,g,c,,a,,
,,b,,d,,f
"""

# Solve the puzzle
print(solve_puzzle(puzzle))
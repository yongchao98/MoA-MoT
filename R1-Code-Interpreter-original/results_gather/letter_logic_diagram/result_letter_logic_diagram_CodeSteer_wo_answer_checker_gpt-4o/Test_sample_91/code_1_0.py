def solve_puzzle(grid):
    from collections import defaultdict

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the letters used in each row and column
    all_letters = set('abcdefg')
    rows = [set(row) - {''} for row in grid]
    cols = [set(grid[r][c] for r in range(7)) - {''} for c in range(7)]

    # Find the minor diagonal indices
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Backtracking function
    def backtrack(r, c):
        if r == 7:  # If we've filled all rows, the puzzle is solved
            return True
        if c == 7:  # Move to the next row
            return backtrack(r + 1, 0)
        if grid[r][c] != '':  # Skip pre-filled cells
            return backtrack(r, c + 1)

        # Determine possible letters for this cell
        possible_letters = all_letters - rows[r] - cols[c]
        if (r, c) in minor_diagonal_indices:
            # If on the minor diagonal, all must be the same letter
            if minor_diagonal_indices.index((r, c)) == 0:
                # Choose a letter for the minor diagonal
                possible_letters &= all_letters - {grid[i][j] for i, j in minor_diagonal_indices if grid[i][j] != ''}
            else:
                # Must match the first minor diagonal letter
                possible_letters &= {grid[minor_diagonal_indices[0][0]][minor_diagonal_indices[0][1]]}

        for letter in possible_letters:
            # Place the letter
            grid[r][c] = letter
            rows[r].add(letter)
            cols[c].add(letter)

            # Recursively attempt to fill the rest of the grid
            if backtrack(r, c + 1):
                return True

            # Backtrack
            grid[r][c] = ''
            rows[r].remove(letter)
            cols[c].remove(letter)

        return False

    # Start backtracking from the top-left corner
    if backtrack(0, 0):
        # Format the output
        result = '\n'.join(','.join(row) for row in grid)
        return f"<<<\n{result}\n>>>"
    else:
        return "No solution found."

# Given puzzle
puzzle = """
d,g,c,e,,a,
g,c,,,,,
,,f,,,d,
e,,,,d,g,
,,,d,g,,e
a,,,,,,f
,,,,e,,a
"""

# Solve the puzzle
print(solve_puzzle(puzzle))
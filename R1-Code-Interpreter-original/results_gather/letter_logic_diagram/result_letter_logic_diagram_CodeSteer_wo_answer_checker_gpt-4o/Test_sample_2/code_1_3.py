def solve_puzzle(grid):
    from collections import defaultdict

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the letters used in each row and column
    all_letters = set('abcdefg')
    rows = [set(row) for row in grid]
    cols = [set(grid[r][c] for r in range(7)) for c in range(7)]

    # Find the minor diagonal indices
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Function to check if placing a letter is valid
    def is_valid(r, c, letter):
        return letter not in rows[r] and letter not in cols[c]

    # Function to print the current state of the grid
    def print_grid():
        for row in grid:
            print(','.join(row))
        print()

    # Backtracking function to fill the grid
    def backtrack(index=0):
        if index == 49:  # All cells are filled
            return True

        r, c = divmod(index, 7)
        if grid[r][c] != '':  # Skip pre-filled cells
            return backtrack(index + 1)

        # Try placing each letter
        for letter in sorted(all_letters, key=lambda l: (l in rows[r], l in cols[c])):
            if is_valid(r, c, letter):
                grid[r][c] = letter
                rows[r].add(letter)
                cols[c].add(letter)

                if backtrack(index + 1):
                    return True

                # Backtrack
                grid[r][c] = ''
                rows[r].remove(letter)
                cols[c].remove(letter)

        return False

    # Try each letter for the minor diagonal
    for diagonal_letter in all_letters:
        if all(is_valid(r, c, diagonal_letter) for r, c in minor_diagonal_indices):
            for r, c in minor_diagonal_indices:
                grid[r][c] = diagonal_letter
                rows[r].add(diagonal_letter)
                cols[c].add(diagonal_letter)

            print(f"Trying diagonal letter: {diagonal_letter}")
            print_grid()

            if backtrack():
                # Format the output
                result = '\n'.join(','.join(row) for row in grid)
                return f"<<<\n{result}\n>>>"

            # Reset the grid for the next attempt
            for r, c in minor_diagonal_indices:
                grid[r][c] = ''
                rows[r].remove(diagonal_letter)
                cols[c].remove(diagonal_letter)

    return "No solution found"

# Given puzzle
puzzle = """
,b,,f,g,,c
b,e,,,a,c,
,f,g,a,c,,b
f,,a,,d,,e
g,,,d,,e,f
a,,,,,f,g
c,d,,e,,g,a
"""

# Solve the puzzle
print(solve_puzzle(puzzle))
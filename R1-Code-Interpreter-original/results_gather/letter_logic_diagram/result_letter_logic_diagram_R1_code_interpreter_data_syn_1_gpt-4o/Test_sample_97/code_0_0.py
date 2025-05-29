def solve_puzzle(grid):
    from itertools import permutations

    # Convert grid to a list of lists for easier manipulation
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Identify the minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Determine possible letters for the minor diagonal
    all_letters = set('abcdefg')
    for letter in all_letters:
        # Check if the letter can be placed on the minor diagonal
        if all(grid[r][c] in ('', letter) for r, c in minor_diagonal_positions):
            # Try to fill the grid with this letter on the minor diagonal
            for r, c in minor_diagonal_positions:
                grid[r][c] = letter

            # Check if the grid can be completed
            if fill_grid(grid):
                return grid

            # Reset the minor diagonal
            for r, c in minor_diagonal_positions:
                grid[r][c] = ''

    return None

def fill_grid(grid):
    # Try to fill the grid row by row
    all_letters = set('abcdefg')
    for r in range(7):
        missing_letters = all_letters - set(grid[r])
        if not missing_letters:
            continue

        # Find empty positions in the row
        empty_positions = [c for c in range(7) if grid[r][c] == '']

        # Try all permutations of the missing letters
        for perm in permutations(missing_letters):
            # Check if the permutation can be placed in the row
            if all(grid[r][c] == '' or grid[r][c] == perm[i] for i, c in enumerate(empty_positions)):
                # Place the permutation
                for i, c in enumerate(empty_positions):
                    grid[r][c] = perm[i]

                # Check if the grid is valid so far
                if all(len(set(grid[r])) == 7 for r in range(7)) and all(len(set(grid[c][r] for c in range(7))) == 7 for r in range(7)):
                    return True

                # Reset the row
                for i, c in enumerate(empty_positions):
                    grid[r][c] = ''

    return False

# Initial grid
grid = """
,a,e,g,f,b,
a,,g,f,,,
,,,,c,d,a
,f,b,c,,a,
,b,c,,a,,g
b,,,,,g,
c,,,e,g,,
"""

# Solve the puzzle
solution = solve_puzzle(grid)

# Format the solution for output
if solution:
    for row in solution:
        print(','.join(row))
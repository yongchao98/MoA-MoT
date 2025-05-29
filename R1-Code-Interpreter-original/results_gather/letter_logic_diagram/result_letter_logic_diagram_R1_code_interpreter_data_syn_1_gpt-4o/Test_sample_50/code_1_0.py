def solve_puzzle(grid):
    from itertools import permutations

    # Convert grid to a list of lists for easier manipulation
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Identify the minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Determine possible letters for the minor diagonal
    all_letters = set('abcdefg')
    possible_letters = all_letters.copy()

    for r, c in minor_diagonal_positions:
        if grid[r][c] != '':
            possible_letters.intersection_update(grid[r][c])

    # Try each possible letter for the minor diagonal
    for letter in possible_letters:
        # Create a copy of the grid to work with
        temp_grid = [row[:] for row in grid]

        # Fill the minor diagonal with the chosen letter
        for r, c in minor_diagonal_positions:
            temp_grid[r][c] = letter

        # Check if the current grid is valid and fill the rest
        if fill_grid(temp_grid):
            return temp_grid

    return None

def fill_grid(grid):
    from itertools import permutations

    # Try to fill the grid row by row
    all_letters = set('abcdefg')
    for r in range(7):
        row_letters = set(grid[r])
        missing_letters = all_letters - row_letters

        # Find empty positions in the row
        empty_positions = [c for c in range(7) if grid[r][c] == '']

        # Try all permutations of missing letters in empty positions
        for perm in permutations(missing_letters):
            for idx, c in enumerate(empty_positions):
                grid[r][c] = perm[idx]

            # Check if the current row is valid
            if all(len(set(grid[r])) == 7 for r in range(7)) and all(len(set(grid[c][r] for c in range(7))) == 7 for r in range(7)):
                break
            else:
                # Reset the row if not valid
                for c in empty_positions:
                    grid[r][c] = ''

    # Check if the entire grid is valid
    return all(len(set(grid[r])) == 7 for r in range(7)) and all(len(set(grid[c][r] for c in range(7))) == 7 for r in range(7))

# Initial grid
grid = """
c,,,,b,,
,f,e,b,,g,
,,b,,,,
,b,a,,c,d,f
b,a,,,,,
,g,,,f,e,
,c,,f,,,
"""

# Solve the puzzle
solution = solve_puzzle(grid)

# Format the solution for output
if solution:
    for row in solution:
        print(','.join(row))
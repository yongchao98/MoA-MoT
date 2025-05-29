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

    # Extract the minor diagonal
    minor_diagonal = [grid[i][6-i] for i in range(7) if grid[i][6-i] != '']

    # If the minor diagonal is not consistent, return None
    if len(set(minor_diagonal)) > 1:
        return None

    # Determine the letter for the minor diagonal
    minor_letter = minor_diagonal[0] if minor_diagonal else None

    # Fill the grid
    def fill_grid(grid, row=0):
        if row == 7:
            return is_valid(grid)

        # Find missing letters for the current row
        missing_letters = set('abcdefg') - set(grid[row])

        # Try all permutations of missing letters
        for perm in permutations(missing_letters):
            new_grid = [list(r) for r in grid]
            idx = 0
            for col in range(7):
                if new_grid[row][col] == '':
                    new_grid[row][col] = perm[idx]
                    idx += 1

            # Check if the minor diagonal condition is satisfied
            if minor_letter and new_grid[row][6-row] != minor_letter:
                continue

            # Recurse to the next row
            if fill_grid(new_grid, row + 1):
                for r in range(7):
                    grid[r] = new_grid[r]
                return True

        return False

    # Start filling the grid
    if fill_grid(grid):
        return grid
    else:
        return None

# Initial grid
grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")
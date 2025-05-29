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

    # Helper function to fill the grid
    def fill_grid(grid):
        # Find the missing letters for each row and column
        rows_missing = [set('abcdefg') - set(row) for row in grid]
        cols_missing = [set('abcdefg') - set(grid[i][j] for i in range(7)) for j in range(7)]

        # Try to fill the grid
        for perm in permutations('abcdefg'):
            # Fill the minor diagonal with the same letter
            for i in range(7):
                grid[i][6-i] = perm[0]

            # Fill the rest of the grid
            for i in range(7):
                for j in range(7):
                    if grid[i][j] == '':
                        possible_letters = rows_missing[i] & cols_missing[j]
                        if len(possible_letters) == 1:
                            grid[i][j] = possible_letters.pop()

            if is_valid(grid):
                return grid

        return None

    # Convert the input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Fill the grid
    filled_grid = fill_grid(grid)

    # Convert the grid back to the required format
    if filled_grid:
        return '\n'.join(','.join(row) for row in filled_grid)
    else:
        return "No solution found"

# Input grid
input_grid = """
,b,,f,g,,c
b,e,,,a,c,
,f,g,a,c,,b
f,,a,,d,,e
g,,,d,,e,f
a,,,,,f,g
c,d,,e,,g,a
"""

# Solve the puzzle
solution = solve_puzzle(input_grid)
print(solution)
def solve_puzzle(grid):
    from itertools import permutations

    # Check if a grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Fill the grid
    def fill_grid(grid):
        # Find the missing letters for each row and column
        all_letters = set('abcdefg')
        rows_missing = [all_letters - set(row) for row in grid]
        cols_missing = [all_letters - set(grid[i][j] for i in range(7)) for j in range(7)]

        # Try to fill the grid
        for perm in permutations('abcdefg'):
            # Set the minor diagonal
            for i in range(7):
                grid[i][6-i] = perm[i]

            # Fill the rest of the grid
            for i in range(7):
                for j in range(7):
                    if grid[i][j] == '':
                        possible_letters = rows_missing[i] & cols_missing[j]
                        if len(possible_letters) == 1:
                            letter = possible_letters.pop()
                            grid[i][j] = letter
                            rows_missing[i].remove(letter)
                            cols_missing[j].remove(letter)

            if is_valid(grid):
                return grid

        return None

    # Convert input grid to a list of lists
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
d,,,f,a,,
b,,f,,,e,d
,f,a,c,e,d,b
f,a,c,,,,
a,,e,,,g,f
c,e,,b,,,a
e,d,,g,f,a,
"""

# Solve the puzzle
solution = solve_puzzle(input_grid)
print(solution)
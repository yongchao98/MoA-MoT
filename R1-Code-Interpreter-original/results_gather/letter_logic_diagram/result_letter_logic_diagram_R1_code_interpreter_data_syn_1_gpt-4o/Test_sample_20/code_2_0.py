def solve_puzzle(grid):
    from itertools import permutations

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Identify the minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Determine the possible letters for the minor diagonal
    all_letters = set('abcdefg')
    used_letters = set(grid[i][j] for i, j in minor_diagonal_positions if grid[i][j])
    possible_diagonal_letters = all_letters - used_letters

    # Try each possible letter for the minor diagonal
    for diagonal_letter in possible_diagonal_letters:
        # Create a copy of the grid to work with
        temp_grid = [row[:] for row in grid]

        # Fill the minor diagonal with the chosen letter
        for i, j in minor_diagonal_positions:
            temp_grid[i][j] = diagonal_letter

        # Try to fill the rest of the grid
        def is_valid():
            for i in range(7):
                if len(set(temp_grid[i])) != 7 or len(set(row[i] for row in temp_grid)) != 7:
                    return False
            return True

        def fill_grid():
            for i in range(7):
                missing_letters = list(all_letters - set(temp_grid[i]))
                for perm in permutations(missing_letters):
                    temp_row = []
                    perm_index = 0
                    for j in range(7):
                        if temp_grid[i][j]:
                            temp_row.append(temp_grid[i][j])
                        else:
                            temp_row.append(perm[perm_index])
                            perm_index += 1
                    if len(set(temp_row)) == 7:
                        temp_grid[i] = temp_row
                        break

        fill_grid()

        if is_valid():
            return temp_grid

    return None

# Given puzzle grid
puzzle_grid = """
b,,,g,,c,e
,d,g,a,,e,
d,,a,c,,b,f
g,,,e,,,d
,,e,b,,,
,e,,f,d,g,a
e,b,f,d,g,,c
"""

# Solve the puzzle
solution = solve_puzzle(puzzle_grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found.")
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

        # Backtracking function to fill the grid
        def is_valid():
            for i in range(7):
                if len(set(temp_grid[i])) != 7 or len(set(row[i] for row in temp_grid)) != 7:
                    return False
            return True

        def backtrack(row, col):
            if row == 7:
                return is_valid()

            if col == 7:
                return backtrack(row + 1, 0)

            if temp_grid[row][col]:
                return backtrack(row, col + 1)

            for letter in all_letters:
                if letter not in temp_grid[row] and letter not in [temp_grid[r][col] for r in range(7)]:
                    temp_grid[row][col] = letter
                    if backtrack(row, col + 1):
                        return True
                    temp_grid[row][col] = ''

            return False

        if backtrack(0, 0):
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
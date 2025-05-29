def solve_puzzle(grid):
    from itertools import permutations

    # Convert the grid into a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(i, 6 - i) for i in range(7)]
    diagonal_letter = None
    for i, j in minor_diagonal_indices:
        if grid[i][j]:
            diagonal_letter = grid[i][j]
            break

    # If no letter is found on the diagonal, choose one arbitrarily (e.g., 'a')
    if not diagonal_letter:
        diagonal_letter = 'a'

    # Fill the minor diagonal with the determined letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = diagonal_letter

    # Function to check if the grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
                return False
        return True

    # Try to fill the grid
    letters = set('abcdefg')
    for perm in permutations(letters):
        for i in range(7):
            missing_letters = letters - set(grid[i])
            if len(missing_letters) == 1:
                grid[i] = [c if c else missing_letters.pop() for c in grid[i]]

        if is_valid(grid):
            break

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Given puzzle
puzzle = """
,,,e,d,c,
f,,e,,,,
,,d,,g,,
e,,c,,b,f,
d,,,b,,a,e
c,g,b,,a,e,d
,b,f,,,d,
"""

# Solve the puzzle
print(solve_puzzle(puzzle))
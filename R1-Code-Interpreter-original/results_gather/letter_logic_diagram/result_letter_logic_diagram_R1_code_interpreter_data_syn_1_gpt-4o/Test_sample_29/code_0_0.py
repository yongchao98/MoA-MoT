def solve_puzzle(grid):
    from collections import Counter

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the missing letters for each row and column
    all_letters = set('abcdefg')
    missing_in_rows = [all_letters - set(row) for row in grid]
    missing_in_cols = [all_letters - set(col) for col in zip(*grid)]

    # Determine the letter for the minor diagonal
    # The minor diagonal indices are (0,6), (1,5), (2,4), (3,3), (4,2), (5,1), (6,0)
    minor_diagonal_indices = [(i, 6-i) for i in range(7)]
    minor_diagonal_letters = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    if minor_diagonal_letters:
        minor_diagonal_letter = Counter(minor_diagonal_letters).most_common(1)[0][0]
    else:
        # If no letters are pre-filled on the diagonal, choose one that can fit
        possible_letters = all_letters - set(minor_diagonal_letters)
        minor_diagonal_letter = possible_letters.pop()

    # Fill in the grid
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                if (i, j) in minor_diagonal_indices:
                    grid[i][j] = minor_diagonal_letter
                else:
                    # Find a letter that can fit in this position
                    possible_letters = missing_in_rows[i] & missing_in_cols[j]
                    grid[i][j] = possible_letters.pop()
                    missing_in_rows[i].remove(grid[i][j])
                    missing_in_cols[j].remove(grid[i][j])

    # Convert grid back to the required format
    return '\n'.join([','.join(row) for row in grid])

# Given puzzle
puzzle = """
a,e,f,g,,c,d
e,f,g,b,,d,a
,g,b,,,,e
,b,,,a,,f
,c,d,,e,,g
c,d,a,e,f,g,b
d,a,e,f,,b,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(solution)
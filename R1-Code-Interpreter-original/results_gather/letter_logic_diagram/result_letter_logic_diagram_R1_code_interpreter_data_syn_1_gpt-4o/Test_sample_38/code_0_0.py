def fill_grid(grid):
    letters = {'a', 'b', 'c', 'd', 'e', 'f', 'g'}
    n = 7

    # Determine the letter for the minor diagonal
    # The minor diagonal indices are (0,6), (1,5), (2,4), (3,3), (4,2), (5,1), (6,0)
    minor_diagonal_letter = None
    for i in range(n):
        if grid[i][n-1-i] != ' ':
            minor_diagonal_letter = grid[i][n-1-i]
            break

    # If no letter is found on the minor diagonal, choose one that is not used
    if not minor_diagonal_letter:
        used_letters = set()
        for i in range(n):
            for j in range(n):
                if grid[i][j] != ' ':
                    used_letters.add(grid[i][j])
        minor_diagonal_letter = (letters - used_letters).pop()

    # Fill the grid
    for i in range(n):
        row_letters = set(grid[i])
        missing_row_letters = letters - row_letters

        for j in range(n):
            if grid[i][j] == ' ':
                if j == n-1-i:  # Minor diagonal
                    grid[i][j] = minor_diagonal_letter
                else:
                    # Find a letter that is not in the current column
                    column_letters = {grid[k][j] for k in range(n)}
                    possible_letters = missing_row_letters - column_letters
                    grid[i][j] = possible_letters.pop()
                    missing_row_letters.remove(grid[i][j])

    return grid

# Initial grid
grid = [
    ['a', 'e', ' ', 'f', 'g', 'c', ' '],
    [' ', 'b', ' ', ' ', ' ', 'd', 'a'],
    ['b', 'f', 'g', 'c', 'd', 'a', 'e'],
    ['f', ' ', 'c', ' ', 'a', ' ', 'b'],
    ['g', 'c', 'd', 'a', 'e', 'b', 'f'],
    ['c', ' ', ' ', 'e', 'b', ' ', 'g'],
    ['d', ' ', ' ', ' ', 'f', 'g', 'c']
]

filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
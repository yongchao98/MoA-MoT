def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        ['b', None, None, 'c', 'g', 'f', 'a'],
        [None, 'e', 'c', 'g', 'f', 'a', 'b'],
        ['e', 'c', 'g', 'f', 'a', 'b', 'd'],
        ['c', None, 'f', None, 'b', None, 'e'],
        ['g', 'f', None, 'b', 'd', None, None],
        [None, None, 'b', None, 'e', 'c', 'g'],
        [None, None, 'd', 'e', None, 'g', 'f']
    ]

    # Determine possible letters for the minor diagonal
    def find_diagonal_letter():
        for letter in 'abcdefg':
            if all((letter not in grid[i] and
                    all(grid[k][6-i] != letter for k in range(7)))
                   for i in range(7)):
                return letter
        return None

    # Choose a letter for the diagonal
    diagonal_letter = find_diagonal_letter()
    if diagonal_letter is None:
        raise ValueError("No valid letter found for the diagonal")

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Helper function to check if a letter can be placed at grid[i][j]
    def can_place(i, j, letter):
        return (letter not in grid[i] and
                all(grid[k][j] != letter for k in range(7)))

    # Backtracking function to fill the grid
    def fill_grid(i=0, j=0):
        if i == 7:
            return True
        if j == 7:
            return fill_grid(i + 1, 0)
        if grid[i][j] is not None:
            return fill_grid(i, j + 1)

        for letter in 'abcdefg':
            if can_place(i, j, letter):
                grid[i][j] = letter
                if fill_grid(i, j + 1):
                    return True
                grid[i][j] = None

        return False

    # Fill the grid
    if not fill_grid():
        raise ValueError("No solution found")

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()
def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid, diagonal_letter):
                            return True
                        grid[row][col] = ''
                return False
    return True

def fill_grid():
    grid = [
        ['f', 'e', 'g', 'a', '', '', ''],
        ['e', 'g', 'a', '', '', '', 'f'],
        ['g', '', '', '', 'b', 'f', 'e'],
        ['', 'd', '', '', 'f', 'e', 'g'],
        ['d', '', 'b', '', 'e', '', 'a'],
        ['', 'b', 'f', 'e', '', 'a', 'd'],
        ['', '', 'e', '', 'a', 'd', 'c']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    for i, j in diagonal_positions:
        if grid[i][j] in letters:
            letters.remove(grid[i][j])
    
    # Choose the letter for the diagonal
    diagonal_letter = letters.pop()

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Solve the grid
    solve(grid, diagonal_letter)

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()
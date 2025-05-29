def fill_grid(grid):
    letters = set('abcdefg')
    n = len(grid)
    
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in letters:
        if all(grid[i][n-1-i] in ('', letter) for i in range(n)):
            minor_diagonal_letter = letter
            break
    
    # Fill the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                if j == n-1-i:  # Minor diagonal
                    grid[i][j] = minor_diagonal_letter
                else:
                    # Determine possible letters for this cell
                    row_letters = set(grid[i])
                    col_letters = set(grid[k][j] for k in range(n))
                    possible_letters = letters - row_letters - col_letters
                    grid[i][j] = possible_letters.pop()
    
    return grid

# Initial grid
grid = [
    ['e', 'f', '', 'd', '', 'a', 'c'],
    ['', 'b', '', '', 'a', '', 'e'],
    ['', 'd', 'g', 'a', '', 'e', ''],
    ['d', 'g', 'a', 'c', 'e', 'f', 'b'],
    ['g', 'a', '', 'e', 'f', 'b', 'd'],
    ['', 'c', '', '', '', 'd', 'g'],
    ['', 'e', 'f', 'b', 'd', 'g', 'a']
]

filled_grid = fill_grid(grid)
for row in filled_grid:
    print(','.join(row))
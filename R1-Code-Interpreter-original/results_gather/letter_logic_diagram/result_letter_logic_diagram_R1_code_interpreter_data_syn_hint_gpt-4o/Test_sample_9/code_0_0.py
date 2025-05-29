def fill_grid(grid):
    letters = set('abcdefg')
    n = len(grid)
    
    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    diagonal_letter = None
    for letter in letters:
        if all(grid[i][n-i-1] in ('', letter) for i in range(n)):
            diagonal_letter = letter
            break
    
    # Fill the minor diagonal
    for i in range(n):
        grid[i][n-i-1] = diagonal_letter
    
    # Fill the rest of the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                possible_letters = letters - row_letters - col_letters
                grid[i][j] = possible_letters.pop()
    
    return grid

# Initial grid setup
grid = [
    ['', 'c', '', '', 'b', '', ''],
    ['c', 'f', '', '', 'a', '', ''],
    ['', '', 'b', 'a', 'e', '', 'c'],
    ['', 'b', '', '', '', '', 'f'],
    ['', '', 'e', 'g', 'c', 'f', ''],
    ['', '', 'g', '', '', 'd', ''],
    ['e', '', '', '', '', 'b', '']
]

filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
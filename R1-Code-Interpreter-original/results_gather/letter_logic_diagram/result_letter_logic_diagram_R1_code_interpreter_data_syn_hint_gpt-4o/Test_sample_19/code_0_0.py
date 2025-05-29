def fill_grid(grid):
    # Determine the letter for the minor diagonal
    letters = set('abcdefg')
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    
    # Check which letter can be placed on the minor diagonal
    for letter in letters:
        can_place = True
        for r, c in minor_diagonal_positions:
            if letter in grid[r] or any(grid[i][c] == letter for i in range(7)):
                can_place = False
                break
        if can_place:
            minor_diagonal_letter = letter
            break
    
    # Place the letter on the minor diagonal
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        missing_letters = letters - set(grid[r])
        for c in range(7):
            if grid[r][c] == '':
                for letter in missing_letters:
                    if all(grid[i][c] != letter for i in range(7)):
                        grid[r][c] = letter
                        missing_letters.remove(letter)
                        break
    
    return grid

# Initial grid setup
grid = [
    ['g', 'b', 'f', '', 'c', 'a', ''],
    ['', 'f', 'd', '', 'a', 'e', ''],
    ['', 'd', 'c', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'e', '', 'b', '', ''],
    ['a', 'e', '', '', '', 'd', 'c'],
    ['e', '', '', '', '', 'c', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
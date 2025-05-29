def fill_grid(grid):
    # Determine the letter for the minor diagonal
    diagonal_letter = 'c'
    
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter
    
    # Function to find missing letters in a row
    def find_missing_letters(row):
        all_letters = set('abcdefg')
        return list(all_letters - set(row))
    
    # Fill the grid row by row
    for i in range(7):
        missing_letters = find_missing_letters(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that is not in the current column
                for letter in missing_letters:
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break
    
    return grid

# Initial grid
grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
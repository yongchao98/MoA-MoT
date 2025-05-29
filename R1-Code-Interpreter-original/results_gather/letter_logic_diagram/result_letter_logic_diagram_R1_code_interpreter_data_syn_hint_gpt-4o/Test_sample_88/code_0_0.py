# Initial grid setup
grid = [
    ['f', '', '', '', 'e', 'a', 'g'],
    ['c', 'd', 'b', '', 'a', 'g', 'f'],
    ['', 'b', '', 'a', 'g', 'f', 'c'],
    ['b', 'e', '', 'g', 'f', 'c', 'd'],
    ['', 'a', 'g', 'f', 'c', 'd', 'b'],
    ['', 'g', 'f', '', '', 'b', 'e'],
    ['g', '', '', '', '', 'e', 'a']
]

# Function to fill the grid
def fill_grid(grid):
    # Determine the letter for the minor diagonal
    # We choose 'f' as it appears in the diagonal already
    diagonal_letter = 'f'
    
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
                grid[i][j] = missing_letters.pop(0)
    
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
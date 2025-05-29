from itertools import permutations

# Initial grid with pre-filled values
grid = [
    ['', '', 'c', 'f', '', '', ''],
    ['e', 'c', 'f', '', '', '', ''],
    ['', 'f', 'g', 'b', 'a', '', 'e'],
    ['f', 'g', 'b', 'a', '', '', 'c'],
    ['g', 'b', 'a', 'd', 'e', '', 'f'],
    ['', '', '', '', 'c', '', ''],
    ['a', '', 'e', '', '', '', '']
]

# Function to check if a letter can be placed on the minor diagonal
def can_place_on_diagonal(letter):
    for i in range(7):
        if grid[i][6-i] not in ('', letter):
            return False
    return True

# Function to fill the grid
def fill_grid(diagonal_letter):
    # Place the diagonal letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter
    
    # Fill each row
    for i in range(7):
        missing_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                for letter in missing_letters:
                    # Check if the letter can be placed in the column
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break

# Try each letter for the diagonal
for letter in 'abcdefg':
    if can_place_on_diagonal(letter):
        fill_grid(letter)
        break

# Print the filled grid
for row in grid:
    print(','.join(row))
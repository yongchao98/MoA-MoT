from itertools import permutations

# Initial grid setup
grid = [
    ['', 'b', 'a', '', '', '', ''],
    ['b', '', '', 'd', '', '', ''],
    ['', '', 'd', 'f', '', '', 'b'],
    ['g', 'd', '', '', 'e', '', ''],
    ['d', 'f', 'c', 'e', 'b', 'a', ''],
    ['', '', '', '', '', '', 'd'],
    ['', 'e', '', '', '', '', 'f']
]

# Function to check if a letter can be placed at a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to fill the grid
def fill_grid(grid):
    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
            
            # Fill the rest of the grid
            for row in range(7):
                for col in range(7):
                    if grid[row][col] == '':
                        for l in 'abcdefg':
                            if can_place(grid, row, col, l):
                                grid[row][col] = l
                                break
            return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))
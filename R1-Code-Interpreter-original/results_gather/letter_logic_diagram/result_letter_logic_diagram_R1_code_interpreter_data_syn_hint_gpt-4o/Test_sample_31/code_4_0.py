from itertools import permutations

def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter):
    # Find the next empty cell
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                # If it's a diagonal cell, it must be the diagonal_letter
                if row + col == 6:
                    if is_valid(grid, row, col, diagonal_letter):
                        grid[row][col] = diagonal_letter
                        if solve(grid, diagonal_letter):
                            return True
                        grid[row][col] = ''
                else:
                    for letter in 'abcdefg':
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if solve(grid, diagonal_letter):
                                return True
                            grid[row][col] = ''
                return False
    return True

# Initial grid with given letters
grid = [
    ['e', 'a', '', '', '', 'c', 'f'],
    ['', 'b', '', '', 'c', 'f', ''],
    ['b', 'g', 'd', 'c', 'f', '', 'a'],
    ['g', 'd', '', 'f', 'e', '', 'b'],
    ['d', '', '', 'e', '', '', 'g'],
    ['', 'f', 'e', '', 'b', '', 'd'],
    ['', '', 'a', '', 'g', 'd', 'c']
]

# Try each letter as the diagonal letter
for diagonal_letter in 'abcdefg':
    if solve(grid, diagonal_letter):
        break

# Print the solved grid
for row in grid:
    print(','.join(row))
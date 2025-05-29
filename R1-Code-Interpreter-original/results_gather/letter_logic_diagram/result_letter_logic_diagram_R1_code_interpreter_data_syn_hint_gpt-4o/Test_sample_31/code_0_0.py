def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    # Check if the letter is on the minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    return True

def solve(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid):
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

# Solve the puzzle
solve(grid)

# Print the solved grid
for row in grid:
    print(','.join(row))
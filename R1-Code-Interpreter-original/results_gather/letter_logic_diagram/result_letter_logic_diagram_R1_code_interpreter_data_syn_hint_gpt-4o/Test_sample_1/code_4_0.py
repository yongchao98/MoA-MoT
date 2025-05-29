def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] is None:
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid):
                            return True
                        grid[row][col] = None
                return False
    return True

# Define the initial grid with empty spaces as None
grid = [
    ['b', None, None, 'c', 'g', 'f', 'a'],
    [None, 'e', 'c', 'g', 'f', 'a', 'b'],
    ['e', 'c', 'g', 'f', 'a', 'b', 'd'],
    ['c', None, 'f', None, 'b', None, 'e'],
    ['g', 'f', None, 'b', 'd', None, None],
    [None, None, 'b', None, 'e', 'c', 'g'],
    [None, None, 'd', 'e', None, 'g', 'f']
]

# Fill the minor diagonal with 'a'
for i in range(7):
    grid[i][6-i] = 'a'

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))
def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the row
    if letter in grid[row]:
        return False
    # Check if the letter can be placed in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] is None:
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[row][col] = None
                return False
    return True

def fill_minor_diagonal(grid):
    for letter in 'abcdefg':
        if all(grid[i][6-i] is None or grid[i][6-i] == letter for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            return letter
    return None

# Initial grid with empty cells as None
grid = [
    [None, None, 'a', None, 'c', 'f', 'b'],
    ['g', 'a', None, None, None, 'b', 'e'],
    [None, None, 'c', 'f', 'b', 'e', 'g'],
    ['d', 'c', 'f', None, 'e', 'g', 'a'],
    [None, 'f', 'b', 'e', 'g', None, None],
    ['f', 'b', 'e', 'g', 'a', 'd', None],
    ['b', 'e', None, 'a', 'd', None, 'f']
]

# Fill the minor diagonal
minor_diagonal_letter = fill_minor_diagonal(grid)

# Solve the grid
if minor_diagonal_letter and solve(grid, minor_diagonal_letter):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")
def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

def fill_minor_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

def find_minor_diagonal_letter(grid):
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

# Initial grid
grid = [
    ['g', 'b', '', '', '', 'a', 'c'],
    ['b', 'e', 'f', 'd', 'a', 'c', 'g'],
    ['', '', '', '', 'c', 'g', ''],
    ['f', '', '', '', '', '', 'e'],
    ['', 'a', 'c', '', 'b', 'e', 'f'],
    ['a', '', 'g', 'b', 'e', 'f', ''],
    ['c', 'g', 'b', '', '', 'd', '']
]

# Determine the letter for the minor diagonal
minor_diagonal_letter = find_minor_diagonal_letter(grid)
if minor_diagonal_letter:
    fill_minor_diagonal(grid, minor_diagonal_letter)

# Solve the puzzle
solve(grid)

# Format the output
result = '\n'.join([','.join(row) for row in grid])
print(result)
def is_valid(grid, row, col, char):
    # Check if char is not in the current row
    if char in grid[row]:
        return False
    # Check if char is not in the current column
    for r in range(7):
        if grid[r][col] == char:
            return False
    return True

def solve(grid, minor_diagonal_char):
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '_':
            return backtrack(row, col + 1)
        
        for char in 'abcdefg':
            if (row + col == 6 and char != minor_diagonal_char) or not is_valid(grid, row, col, char):
                continue
            grid[row][col] = char
            if backtrack(row, col + 1):
                return True
            grid[row][col] = '_'
        return False

    return backtrack(0, 0)

# Initial grid with underscores for empty cells
grid = [
    ['c', 'g', '_', 'd', 'e', 'b', '_'],
    ['g', 'a', 'd', 'e', 'b', '_', 'c'],
    ['_', 'd', 'e', '_', 'f', '_', 'g'],
    ['d', 'e', '_', 'f', '_', '_', '_'],
    ['e', '_', '_', 'c', 'g', 'a', 'd'],
    ['b', 'f', 'c', 'g', 'a', '_', '_'],
    ['f', 'c', 'g', 'a', '_', 'e', 'b']
]

# Try each letter for the minor diagonal
for char in 'abcdefg':
    if all(grid[i][6-i] in ('_', char) for i in range(7)):
        if solve(grid, char):
            break

# Print the completed grid
for row in grid:
    print(','.join(row))
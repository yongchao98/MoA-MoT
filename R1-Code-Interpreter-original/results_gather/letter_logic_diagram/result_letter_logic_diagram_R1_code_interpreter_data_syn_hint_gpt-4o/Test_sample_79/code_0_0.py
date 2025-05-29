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
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for char in 'abcdefg':
                    if is_valid(grid, row, col, char):
                        # Check minor diagonal constraint
                        if row + col == 6 and char != minor_diagonal_char:
                            continue
                        grid[row][col] = char
                        if solve(grid, minor_diagonal_char):
                            return True
                        grid[row][col] = ''
                return False
    return True

# Initial grid setup
grid = [
    ['d', 'e', 'c', 'a', 'f', '', 'g'],
    ['', 'c', '', '', 'b', 'g', ''],
    ['c', 'a', 'f', '', '', 'd', 'e'],
    ['a', '', '', 'g', 'd', 'e', 'c'],
    ['f', '', 'g', 'd', 'e', 'c', 'a'],
    ['b', 'g', '', '', 'c', '', 'f'],
    ['g', 'd', 'e', 'c', '', 'f', '']
]

# Determine the minor diagonal character
minor_diagonal_char = grid[0][6]  # The top-right corner

# Solve the puzzle
solve(grid, minor_diagonal_char)

# Print the solution
for row in grid:
    print(','.join(row))
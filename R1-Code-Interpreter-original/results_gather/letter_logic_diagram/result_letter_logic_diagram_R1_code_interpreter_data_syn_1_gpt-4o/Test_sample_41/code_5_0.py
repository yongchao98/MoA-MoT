def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                if (row + col == 6) and (letter != minor_diagonal_letter):
                    continue
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    backtrack(0, 0)
    return grid

# Initial grid setup
grid = [
    ['', 'a', 'g', '', '', '', ''],
    ['a', 'g', '', 'f', 'b', '', 'd'],
    ['', 'c', '', '', '', '', ''],
    ['', 'f', 'b', '', 'd', '', ''],
    ['', '', '', '', '', 'g', ''],
    ['', 'e', '', '', 'g', 'c', ''],
    ['', '', 'a', 'g', 'c', '', 'b']
]

# Determine the minor diagonal letter
minor_diagonal_letter = 'g'  # Based on the given grid, 'g' is a candidate

# Solve the puzzle
solved_grid = solve(grid, minor_diagonal_letter)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))
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
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        if grid[i][6-i] == '' or grid[i][6-i] == minor_diagonal_letter:
            grid[i][6-i] = minor_diagonal_letter
        else:
            return False  # If the diagonal already has a different letter, return False

    # Start backtracking from the first cell
    return backtrack(0, 0)

# Initial grid
grid = [
    ['a', 'b', '', 'g', 'e', 'c', 'd'],
    ['', 'f', 'g', 'e', '', 'd', ''],
    ['f', 'g', 'e', '', '', 'a', ''],
    ['', 'e', 'c', 'd', '', '', 'f'],
    ['e', '', 'd', '', '', 'f', ''],
    ['c', '', 'a', 'b', '', 'g', ''],
    ['', '', 'b', '', 'g', 'e', 'c']
]

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    # Make a deep copy of the grid
    grid_copy = [row[:] for row in grid]
    if solve(grid_copy, letter):
        for row in grid_copy:
            print(','.join(row))
        break
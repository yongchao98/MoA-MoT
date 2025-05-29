def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check if the letter is valid for the minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Successfully filled the grid
    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row
    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle(grid):
    # Fill the minor diagonal with a common letter
    minor_diagonal_letter = 'g'  # Choose 'g' as it is not in any minor diagonal position
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    if solve(grid):
        result = '\n'.join([','.join(row) for row in grid])
        return f"<<<\n{result}\n>>>"
    else:
        return "No solution found."

# Initial grid with empty cells as ''
grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
]

# Solve the puzzle
print(solve_puzzle(grid))
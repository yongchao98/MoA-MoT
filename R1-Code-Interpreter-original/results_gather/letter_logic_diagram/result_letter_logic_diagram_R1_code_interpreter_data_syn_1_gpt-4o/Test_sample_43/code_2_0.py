def is_valid(grid, r, c, letter):
    # Check if the letter can be placed at grid[r][c]
    for i in range(7):
        if grid[r][i] == letter or grid[i][c] == letter:
            return False
    return True

def solve_puzzle(grid, minor_diagonal_letter, r=0, c=0):
    if r == 7:
        return True  # Puzzle solved

    # Move to the next row if we reach the end of a row
    if c == 7:
        return solve_puzzle(grid, minor_diagonal_letter, r + 1, 0)

    # Skip pre-filled cells
    if grid[r][c] != '':
        return solve_puzzle(grid, minor_diagonal_letter, r, c + 1)

    # Determine if this cell is on the minor diagonal
    if (r, c) in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        if is_valid(grid, r, c, minor_diagonal_letter):
            grid[r][c] = minor_diagonal_letter
            if solve_puzzle(grid, minor_diagonal_letter, r, c + 1):
                return True
            grid[r][c] = ''  # Backtrack
    else:
        # Try placing each letter from 'a' to 'g'
        for letter in 'abcdefg':
            if is_valid(grid, r, c, letter):
                grid[r][c] = letter
                if solve_puzzle(grid, minor_diagonal_letter, r, c + 1):
                    return True
                grid[r][c] = ''  # Backtrack

    return False

# Initial grid with empty cells as ''
grid = [
    ['d', 'a', '', '', 'g', '', ''],
    ['', '', '', 'g', 'e', '', ''],
    ['f', '', 'g', '', 'c', '', ''],
    ['', 'g', '', '', '', 'a', ''],
    ['g', 'e', '', '', '', 'f', ''],
    ['e', '', '', '', '', 'b', ''],
    ['', 'd', '', 'f', 'b', '', 'e']
]

# Determine the letter for the minor diagonal
minor_diagonal_letter = 'c'  # Choose a letter that can be used for the minor diagonal

# Solve the puzzle
if solve_puzzle(grid, minor_diagonal_letter):
    # Print the solved grid
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")
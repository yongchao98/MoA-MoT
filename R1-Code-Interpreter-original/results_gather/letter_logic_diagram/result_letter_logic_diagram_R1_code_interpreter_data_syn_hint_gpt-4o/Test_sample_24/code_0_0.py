from collections import Counter

# Initial grid setup
grid = [
    ['f', 'd', 'a', '', '', 'g', ''],
    ['d', '', '', 'b', '', '', 'f'],
    ['', 'c', 'b', '', 'e', '', ''],
    ['', '', 'g', 'e', '', '', ''],
    ['', 'g', 'e', 'f', '', 'a', 'c'],
    ['', 'e', '', '', '', '', ''],
    ['', '', '', 'a', '', '', '']
]

# Function to check if a letter can be placed in a specific position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Determine the letter for the minor diagonal
def find_minor_diagonal_letter(grid):
    # Count occurrences of each letter in the minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letter_count = Counter(grid[r][c] for r, c in minor_diagonal_positions if grid[r][c])
    # Find a letter that can be placed in all minor diagonal positions
    for letter in 'abcdefg':
        if all(can_place(grid, r, c, letter) for r, c in minor_diagonal_positions):
            return letter
    return None

# Fill the minor diagonal
minor_diagonal_letter = find_minor_diagonal_letter(grid)
minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
for r, c in minor_diagonal_positions:
    grid[r][c] = minor_diagonal_letter

# Fill the rest of the grid
def fill_grid(grid):
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                for letter in 'abcdefg':
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        break

fill_grid(grid)

# Print the completed grid
for row in grid:
    print(','.join(row))
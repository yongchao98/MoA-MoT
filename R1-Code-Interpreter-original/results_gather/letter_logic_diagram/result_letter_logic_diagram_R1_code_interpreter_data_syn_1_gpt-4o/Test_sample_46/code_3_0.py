from collections import Counter

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

# Initial grid with pre-filled values
grid = [
    ['', 'g', 'e', '', 'b', '', ''],
    ['g', 'e', 'f', 'b', '', '', ''],
    ['', '', '', '', '', '', 'g'],
    ['f', '', 'd', '', '', 'g', 'e'],
    ['', 'd', '', '', '', '', ''],
    ['d', 'c', '', '', '', '', 'b'],
    ['c', '', 'g', '', '', 'b', '']
]

# Determine the letter for the minor diagonal
# Count occurrences of each letter
letter_count = Counter()
for row in grid:
    for letter in row:
        if letter:
            letter_count[letter] += 1

# Find the letter that can be used for the minor diagonal
# It should be the one with the least occurrences
all_letters = set('abcdefg')
used_letters = set(letter_count.keys())
available_letters = all_letters - used_letters
if available_letters:
    diagonal_letter = available_letters.pop()
else:
    diagonal_letter = min(letter_count, key=letter_count.get)

# Fill the minor diagonal
for i in range(7):
    grid[i][6-i] = diagonal_letter

# Solve the grid
solve(grid)

# Print the completed grid
for row in grid:
    print(','.join(row))
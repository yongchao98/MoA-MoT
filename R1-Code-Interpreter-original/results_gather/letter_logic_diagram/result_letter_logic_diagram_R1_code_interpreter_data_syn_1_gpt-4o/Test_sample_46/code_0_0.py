from collections import Counter

# Initial grid with pre-filled values
grid = [
    ['a', 'g', 'e', '', 'b', '', ''],
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

# Fill the rest of the grid
for i in range(7):
    row_letters = set(grid[i])
    missing_letters = all_letters - row_letters
    for j in range(7):
        if grid[i][j] == '':
            # Find a letter that can be placed here
            for letter in missing_letters:
                # Check if the letter can be placed in this column
                if all(grid[k][j] != letter for k in range(7)):
                    grid[i][j] = letter
                    missing_letters.remove(letter)
                    break

# Print the completed grid
for row in grid:
    print(','.join(row))
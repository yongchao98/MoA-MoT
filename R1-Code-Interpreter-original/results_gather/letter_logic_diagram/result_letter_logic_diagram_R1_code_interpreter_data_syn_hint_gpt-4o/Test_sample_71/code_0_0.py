# Initial grid with pre-filled letters
grid = [
    ['d', '', '', 'b', 'e', 'g'],
    ['', '', 'b', '', 'g', ''],
    ['f', '', '', 'e', 'g', '', ''],
    ['a', '', '', '', '', 'd', ''],
    ['', '', 'g', '', '', '', 'a'],
    ['', 'g', '', 'd', '', '', ''],
    ['', 'c', 'd', '', 'a', 'b', '']
]

# Count occurrences of each letter
from collections import Counter

# Flatten the grid and count occurrences
flat_grid = [cell for row in grid for cell in row if cell]
letter_count = Counter(flat_grid)

# Determine the letter for the minor diagonal
# Choose the letter with the least occurrences
letters = set('abcdefg')
least_common_letter = min(letters, key=lambda x: letter_count.get(x, 0))

print("Least common letter for the minor diagonal:", least_common_letter)
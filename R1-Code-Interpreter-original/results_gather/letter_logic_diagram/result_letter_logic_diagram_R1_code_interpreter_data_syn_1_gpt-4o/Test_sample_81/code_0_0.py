from collections import Counter

# Initial grid with empty cells as None
grid = [
    [None, 'b', 'd', None, 'e', 'g', 'f'],
    [None, None, 'a', 'e', 'g', None, None],
    ['d', 'a', None, None, None, 'c', None],
    ['a', 'e', None, None, None, None, 'b'],
    [None, None, 'f', 'c', 'b', None, 'a'],
    [None, 'f', 'c', None, 'd', None, None],
    ['f', None, 'b', None, None, None, 'g']
]

# Function to find the missing letter for the minor diagonal
def find_minor_diagonal_letter(grid):
    # Check which letters are already on the diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = [grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None]
    letter_counts = Counter(diagonal_letters)
    
    # Find a letter that can be used for the entire diagonal
    for letter in 'abcdefg':
        if letter_counts[letter] == len(diagonal_letters):
            return letter
    return None

# Find the letter for the minor diagonal
minor_diagonal_letter = find_minor_diagonal_letter(grid)

# Fill the minor diagonal with the found letter
for i, j in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
    grid[i][j] = minor_diagonal_letter

# Function to fill the grid
def fill_grid(grid):
    letters = set('abcdefg')
    for i in range(7):
        # Fill rows
        row_letters = set(grid[i])
        missing_row_letters = letters - row_letters
        for j in range(7):
            if grid[i][j] is None:
                grid[i][j] = missing_row_letters.pop()
        
        # Fill columns
        for j in range(7):
            col_letters = set(grid[k][j] for k in range(7))
            missing_col_letters = letters - col_letters
            for k in range(7):
                if grid[k][j] is None:
                    grid[k][j] = missing_col_letters.pop()

fill_grid(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))
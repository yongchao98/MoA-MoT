# Define the grid size
grid_size = 20

# List of known checkmark positions (row, column)
checkmarks = [
    (1, 8), (2, 13), (3, 16), (4, 17), (5, 18),
    (6, 10), (7, 19), (8, 2), (9, 9), (10, 5),
    (11, 15), (12, 12), (13, 19), (14, 3), (15, 10),
    (16, 4), (17, 5), (18, 7), (19, 13)
]

# Initialize the grid with '❌'
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Place the known checkmarks
for row, col in checkmarks:
    grid[row - 1][col - 1] = '✅'

# Function to find the missing checkmark position
def find_missing_checkmark():
    # Iterate over each row to find the missing checkmark
    for row in range(grid_size):
        if '✅' not in grid[row]:
            return row + 1

# Find the missing row
missing_row = find_missing_checkmark()

# Based on the pattern, determine the column for the missing checkmark
# We will assume the pattern is cyclic and use the known positions to predict the missing one
# Since the pattern is not immediately obvious, we will use a simple heuristic
# Let's assume the missing checkmark is in the column that maintains the pattern
# We will use the difference between consecutive checkmarks to predict the next one

# Calculate the differences between consecutive checkmark columns
differences = [(checkmarks[i][1] - checkmarks[i - 1][1]) % grid_size for i in range(1, len(checkmarks))]

# Find the most common difference
from collections import Counter
most_common_difference = Counter(differences).most_common(1)[0][0]

# Predict the column for the missing checkmark
last_known_col = checkmarks[-1][1]
predicted_col = (last_known_col + most_common_difference) % grid_size

# Adjust for 1-based index
if predicted_col == 0:
    predicted_col = grid_size

# Place the missing checkmark
grid[missing_row - 1][predicted_col - 1] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))
# Define the grid
grid = [
    [-5, 3, 2, 10, -6],
    [4, 7, 7, 2, 4],
    [2, -7, -3, 7, 3],
    [-1, 8, 9, 8, -9],
    [-6, 2, -3, 9, 10]
]

# Initialize the selected positions with the given starting point
selected_positions = [(0, 4)]

# Initialize counters for negative numbers in each row and column
row_neg_count = [0] * 5
col_neg_count = [0] * 5

# Update the counters for the initial selection
row_neg_count[0] += 1
col_neg_count[4] += 1

# Function to check if adding a number keeps the product positive
def can_add(row, col):
    # Check if adding this number will keep the product of the row and column positive
    new_row_neg_count = row_neg_count[row] + (1 if grid[row][col] < 0 else 0)
    new_col_neg_count = col_neg_count[col] + (1 if grid[row][col] < 0 else 0)
    return new_row_neg_count % 2 == 0 and new_col_neg_count % 2 == 0

# Try to select numbers
for r in range(5):
    for c in range(5):
        if len(selected_positions) >= 10:
            break
        if (r, c) not in selected_positions and can_add(r, c):
            selected_positions.append((r, c))
            if grid[r][c] < 0:
                row_neg_count[r] += 1
                col_neg_count[c] += 1

# Print the selected positions
print(selected_positions)
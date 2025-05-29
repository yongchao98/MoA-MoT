# Define the grid
grid = [
    [4, 7, -7, 6, 2],
    [-6, -3, 3, 10, 2],
    [9, 4, 10, 5, 4],
    [-4, 4, 8, 5, 7],
    [8, 1, 10, -2, 10]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 4)]

# Initialize negative counts for rows and columns
row_neg_count = [0] * 5
col_neg_count = [0] * 5

# Update the negative count for the initial selection
if grid[0][4] < 0:
    row_neg_count[0] += 1
    col_neg_count[4] += 1

# Function to check if adding a number keeps the product positive
def can_add_number(r, c):
    # Temporarily update negative counts
    temp_row_neg = row_neg_count[:]
    temp_col_neg = col_neg_count[:]
    if grid[r][c] < 0:
        temp_row_neg[r] += 1
        temp_col_neg[c] += 1
    # Check if all rows and columns have even negative counts
    for i in range(5):
        if temp_row_neg[i] % 2 != 0 or temp_col_neg[i] % 2 != 0:
            return False
    return True

# Select numbers
for r in range(5):
    for c in range(5):
        if (r, c) != (0, 4) and len(selected_positions) < 10:
            if can_add_number(r, c):
                selected_positions.append((r, c))
                # Update the actual negative counts
                if grid[r][c] < 0:
                    row_neg_count[r] += 1
                    col_neg_count[c] += 1

# Print the selected positions
print(selected_positions)
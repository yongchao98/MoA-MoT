# Define the grid
grid = [
    [-5, 3, 2, 10, -6],
    [4, 7, 7, 2, 4],
    [2, -7, -3, 7, 3],
    [-1, 8, 9, 8, -9],
    [-6, 2, -3, 9, 10]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 4)]

# Function to check if a number can be selected
def can_select(row, col, selected_positions):
    # Check if the position is already selected
    if (row, col) in selected_positions:
        return False
    # Check if selecting this number will keep the product positive
    row_neg_count = sum(1 for r, c in selected_positions if r == row and grid[r][c] < 0)
    col_neg_count = sum(1 for r, c in selected_positions if c == col and grid[r][c] < 0)
    current_neg_count = (grid[row][col] < 0)
    # Ensure even number of negatives in row and column
    if (row_neg_count + current_neg_count) % 2 != 0 or (col_neg_count + current_neg_count) % 2 != 0:
        return False
    return True

# Select numbers from the grid
for row in range(5):
    for col in range(5):
        if len(selected_positions) < 10 and can_select(row, col, selected_positions):
            selected_positions.append((row, col))

# Format the output
output = ', '.join(f"{r} {c}" for r, c in selected_positions)
print(f"<<<{output}>>>")
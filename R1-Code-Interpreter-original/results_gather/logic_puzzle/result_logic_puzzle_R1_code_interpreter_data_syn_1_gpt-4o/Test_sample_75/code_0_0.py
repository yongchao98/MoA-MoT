import itertools

# Define the grid
grid = [
    [2, -6, -10, 1, 2],
    [-3, 2, 7, -3, 6],
    [3, 7, 4, 9, 1],
    [-4, 4, 3, 7, 7],
    [8, 4, 3, 5, 10]
]

# Initial selection
selected_positions = [(0, 4)]

# Function to check if the product condition is satisfied
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    # Check if all rows and columns have an even number of negative numbers
    return all(count % 2 == 0 for count in row_count) and all(count % 2 == 0 for count in col_count)

# Try to select 10 numbers
for r in range(5):
    for c in range(5):
        if (r, c) != (0, 4) and len(selected_positions) < 10:
            selected_positions.append((r, c))
            if not is_valid_selection(selected_positions):
                selected_positions.pop()

# Ensure we have exactly 10 numbers
while len(selected_positions) < 10:
    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if is_valid_selection(selected_positions):
                    break
                else:
                    selected_positions.pop()
        if len(selected_positions) == 10:
            break

# Sort the selected positions
selected_positions.sort()

# Print the selected positions
print(selected_positions)
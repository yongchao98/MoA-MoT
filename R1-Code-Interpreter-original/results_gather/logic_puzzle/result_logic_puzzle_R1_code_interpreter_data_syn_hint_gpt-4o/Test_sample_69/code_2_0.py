# Define the grid
grid = [
    [2, -10, -3, 1, -1],
    [5, 6, 6, -9, 3],
    [9, -4, -7, 10, -5],
    [10, -2, 4, 1, -10],
    [-8, 6, -2, 10, 1]
]

# Initialize the selected positions with the given starting point
selected_positions = [(0, 4)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 5
    col_neg_count = [0] * 5
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    # Check if all rows and columns have an even number of negative numbers
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Function to find a valid selection using a greedy approach
def find_selection(selected_positions):
    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if len(selected_positions) == 10 and is_valid_selection(selected_positions):
                    return selected_positions
                if len(selected_positions) < 10:
                    result = find_selection(selected_positions)
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find the valid selection
result = find_selection(selected_positions)

# Print the selected positions
if result:
    print(', '.join(f"{r} {c}" for r, c in result))
else:
    print("No valid selection found")
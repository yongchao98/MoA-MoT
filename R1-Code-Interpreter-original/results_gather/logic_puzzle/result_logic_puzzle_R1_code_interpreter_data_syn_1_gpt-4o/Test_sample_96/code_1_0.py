# Define the grid
grid = [
    [-9, 8, 3, 13, 5, -1, 6],
    [13, -6, -4, 6, 13, 6, -10],
    [7, -6, -2, 7, 6, -9, 3],
    [6, -5, 1, 4, -10, 2, 9],
    [-9, 14, 12, 1, -10, 9, -2],
    [-3, -2, 14, 14, 2, -8, 5],
    [-2, -1, -8, -3, 7, -6, 9]
]

# Initial selected position
selected_positions = [(0, 6)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_count = [0] * 7
    col_count = [0] * 7
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    return all(count % 2 == 0 for count in row_count) and all(count % 2 == 0 for count in col_count)

# Function to select numbers iteratively
def select_numbers():
    all_positions = [(r, c) for r in range(7) for c in range(7)]
    positive_positions = [(r, c) for r, c in all_positions if grid[r][c] > 0]
    negative_positions = [(r, c) for r, c in all_positions if grid[r][c] < 0]

    # Start with positive numbers
    current_selection = selected_positions[:]
    for pos in positive_positions:
        if pos not in current_selection:
            current_selection.append(pos)
            if len(current_selection) == 16:
                if is_valid_selection(current_selection):
                    return current_selection
                current_selection.pop()

    # Balance with negative numbers
    for pos in negative_positions:
        if pos not in current_selection:
            current_selection.append(pos)
            if len(current_selection) == 16:
                if is_valid_selection(current_selection):
                    return current_selection
                current_selection.pop()

    return None

# Find the valid selection
result = select_numbers()

# Print the result
if result:
    print(", ".join(f"{r} {c}" for r, c in result))
else:
    print("No valid selection found")
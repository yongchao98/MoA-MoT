grid = [
    [8, 12, 9, 3, -11, 9],
    [12, -12, -4, 11, 10, 7],
    [-8, 6, 9, 10, 1, 3],
    [5, 11, 6, 7, -9, 11],
    [6, 1, 11, 12, 4, -3],
    [2, 9, 4, 3, 10, 9]
]

# Start with the given number
selected_positions = [(0, 5)]
selected_numbers = [grid[0][5]]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_products = [1] * 6
    col_products = [1] * 6
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Function to find a valid selection
def find_selection():
    # Start with the given number
    selected_positions = [(0, 5)]
    selected_numbers = [grid[0][5]]
    
    # Try to select numbers strategically
    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                selected_numbers.append(grid[r][c])
                if len(selected_positions) == 12:
                    if is_valid_selection(selected_positions):
                        return selected_positions
                elif is_valid_selection(selected_positions):
                    continue
                selected_positions.pop()
                selected_numbers.pop()
    return None

# Find the selection
result = find_selection()
print(result)
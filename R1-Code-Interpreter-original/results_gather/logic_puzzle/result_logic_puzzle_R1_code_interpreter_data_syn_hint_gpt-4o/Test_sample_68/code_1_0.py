grid = [
    [9, 7, 4, 5, 5],
    [-7, 4, -3, -4, 5],
    [7, 1, 3, 2, 8],
    [4, 9, 2, 1, 6],
    [-4, 6, -9, 7, 3]
]

# Start with the given number
selected_positions = [(0, 4)]
selected_numbers = [grid[0][4]]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    # Check if each row and column has an even number of negative numbers
    for count in row_count + col_count:
        if count % 2 != 0:
            return False
    return True

# Recursive function to select numbers
def select_numbers(selected_positions, selected_numbers):
    if len(selected_positions) == 10:
        if is_valid_selection(selected_positions):
            return selected_positions
        else:
            return None

    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                selected_numbers.append(grid[r][c])
                result = select_numbers(selected_positions, selected_numbers)
                if result is not None:
                    return result
                selected_positions.pop()
                selected_numbers.pop()
    return None

# Find a valid selection
result = select_numbers(selected_positions, selected_numbers)

# Output the selected positions
print(result)
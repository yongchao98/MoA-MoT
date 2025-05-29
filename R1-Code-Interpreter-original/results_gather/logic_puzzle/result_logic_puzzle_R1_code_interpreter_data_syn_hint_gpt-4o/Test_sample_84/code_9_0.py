# Define the grid
grid = [
    [12, 2, 8, 9, 4, 10],
    [-9, -12, 10, 2, -8, 7],
    [-4, 9, 4, -8, -10, -11],
    [-4, 12, -1, 11, -5, -2],
    [-1, 8, 11, 3, -6, 2],
    [11, -6, -12, -11, 4, 4]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    # Check each row
    for i in range(6):
        row_numbers = [grid[i][j] for j in range(6) if (i, j) in selected_positions]
        if row_numbers and not is_positive_product(row_numbers):
            return False
    # Check each column
    for j in range(6):
        col_numbers = [grid[i][j] for i in range(6) if (i, j) in selected_positions]
        if col_numbers and not is_positive_product(col_numbers):
            return False
    return True

# Function to check if the product of a list is positive
def is_positive_product(numbers):
    negative_count = sum(1 for num in numbers if num < 0)
    return negative_count % 2 == 0

# Backtracking function to find a valid selection
def find_selection(selected_positions):
    if len(selected_positions) == 12:
        return selected_positions if is_valid_selection(selected_positions) else None
    
    for i in range(6):
        for j in range(6):
            if (i, j) not in selected_positions:
                selected_positions.append((i, j))
                if is_valid_selection(selected_positions):
                    result = find_selection(selected_positions)
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find the valid selection
result = find_selection(selected_positions)

# Print the selected positions
print(result)
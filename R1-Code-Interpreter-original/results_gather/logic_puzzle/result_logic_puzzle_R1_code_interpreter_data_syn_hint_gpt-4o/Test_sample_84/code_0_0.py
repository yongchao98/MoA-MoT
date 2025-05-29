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

# Function to check if the product of a list is positive
def is_positive_product(numbers):
    product = 1
    for num in numbers:
        product *= num
    return product > 0

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

# Try to select 12 numbers
for i in range(6):
    for j in range(6):
        if (i, j) not in selected_positions:
            selected_positions.append((i, j))
            if len(selected_positions) == 12 and is_valid_selection(selected_positions):
                break
            if len(selected_positions) > 12 or not is_valid_selection(selected_positions):
                selected_positions.pop()
    if len(selected_positions) == 12:
        break

# Print the selected positions
print(selected_positions)
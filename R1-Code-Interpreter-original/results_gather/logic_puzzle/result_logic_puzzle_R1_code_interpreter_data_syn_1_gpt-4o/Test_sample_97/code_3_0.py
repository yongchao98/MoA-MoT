# Define the grid
grid = [
    [10, -10, -5, -13, -7, 5, 7],
    [13, 14, 3, 3, 3, -9, 6],
    [2, -10, -9, 5, -3, 12, -13],
    [-3, 14, 12, 9, 11, 6, 5],
    [10, 1, 13, 9, 8, 6, 3],
    [8, 2, -4, 8, 12, 8, 4],
    [13, 1, 14, -2, 1, -13, 1]
]

# Initial selected position
selected_positions = {(0, 6)}

# Function to check if the product of numbers in a list is positive
def is_product_positive(numbers):
    product = 1
    for num in numbers:
        product *= num
    return product > 0

# Function to select numbers row by row
def select_numbers():
    # Start with the initial position
    selected_positions = {(0, 6)}
    selected_count = 1

    # Iterate over each row
    for row in range(7):
        if selected_count >= 16:
            break

        # Select numbers from the current row
        row_numbers = grid[row]
        row_selected = []

        for col in range(7):
            if selected_count >= 16:
                break

            # Check if adding this number keeps the product positive
            if (row, col) not in selected_positions and is_product_positive(row_selected + [row_numbers[col]]):
                row_selected.append(row_numbers[col])
                selected_positions.add((row, col))
                selected_count += 1

        # Ensure at least one number is selected from each row if needed
        if not row_selected and selected_count < 16:
            for col in range(7):
                if (row, col) not in selected_positions:
                    selected_positions.add((row, col))
                    selected_count += 1
                    break

    return list(selected_positions)[:16]  # Ensure only 16 positions are returned

# Get the selected positions
selected_positions = select_numbers()

# Print the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))
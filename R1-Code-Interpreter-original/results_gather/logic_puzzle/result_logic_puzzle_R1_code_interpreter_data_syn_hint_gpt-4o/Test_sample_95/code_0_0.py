# Define the grid
grid = [
    [2, 8, -12, -1, 9, 7],
    [2, 7, 6, 3, 1, -1],
    [9, 12, 9, 8, -10, 2],
    [6, 10, -5, -12, 5, -7],
    [-12, -8, 3, 8, 12, 4],
    [-5, -9, 11, -12, 5, -8]
]

# Initialize the selected positions with the given starting point
selected_positions = [(0, 5)]

# Function to check if the product of a list of numbers is positive
def is_product_positive(numbers):
    negative_count = sum(1 for num in numbers if num < 0)
    return negative_count % 2 == 0

# Select numbers ensuring the product of each row and column is positive
for i in range(6):
    for j in range(6):
        if (i, j) not in selected_positions:
            # Check if adding this number keeps the product positive
            row_numbers = [grid[i][col] for col in range(6) if (i, col) in selected_positions or col == j]
            col_numbers = [grid[row][j] for row in range(6) if (row, j) in selected_positions or row == i]
            
            if is_product_positive(row_numbers) and is_product_positive(col_numbers):
                selected_positions.append((i, j))
            
            # Stop if we have selected 12 numbers
            if len(selected_positions) == 12:
                break
    if len(selected_positions) == 12:
        break

# Format the output
output = ', '.join(f"{r} {c}" for r, c in selected_positions)
print(f"<<<{output}>>>")
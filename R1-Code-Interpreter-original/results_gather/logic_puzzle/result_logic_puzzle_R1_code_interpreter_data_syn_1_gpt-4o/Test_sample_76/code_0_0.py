grid = [
    [-8, 10, -8, -3, 2],
    [4, 10, 5, 2, 9],
    [-5, -4, 5, 9, -2],
    [1, 1, 5, -4, -9],
    [-9, -5, 1, 10, 3]
]

# Start with the given number
selected_positions = [(0, 4)]
selected_numbers = [2]

# Function to check if the product of numbers in a list is positive
def is_product_positive(numbers):
    negative_count = sum(1 for num in numbers if num < 0)
    return negative_count % 2 == 0

# Try to select numbers ensuring the product of each row and column is positive
for i in range(5):
    for j in range(5):
        if len(selected_positions) == 10:
            break
        if (i, j) not in selected_positions:
            # Check if adding this number keeps the product positive
            row_numbers = [grid[i][col] for col in range(5) if (i, col) in selected_positions or col == j]
            col_numbers = [grid[row][j] for row in range(5) if (row, j) in selected_positions or row == i]
            if is_product_positive(row_numbers) and is_product_positive(col_numbers):
                selected_positions.append((i, j))
                selected_numbers.append(grid[i][j])

# Output the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))
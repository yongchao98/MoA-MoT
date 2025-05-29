grid = [
    [9, 5, 1, 12, -2, 11],
    [-5, 6, 6, 11, 10, 8],
    [8, -1, -6, 1, -8, 5],
    [1, 12, 11, 10, -9, 12],
    [3, 3, -12, -4, 5, 12],
    [-1, -10, -3, -2, 4, 12]
]

# Start with the given number
selected_positions = [(0, 5)]
selected_numbers = [11]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(numbers):
    product = 1
    for num in numbers:
        product *= num
    return product > 0

# Function to count negative numbers in a list
def count_negatives(numbers):
    return sum(1 for num in numbers if num < 0)

# Select numbers ensuring the product of each row and column is positive
for i in range(6):
    for j in range(6):
        if len(selected_positions) >= 12:
            break
        if (i, j) not in selected_positions:
            # Check if adding this number keeps the product positive
            row_numbers = [grid[i][col] for row, col in selected_positions if row == i] + [grid[i][j]]
            col_numbers = [grid[row][j] for row, col in selected_positions if col == j] + [grid[i][j]]
            if (count_negatives(row_numbers) % 2 == 0) and (count_negatives(col_numbers) % 2 == 0):
                selected_positions.append((i, j))
                selected_numbers.append(grid[i][j])

# Output the selected positions
print(selected_positions)
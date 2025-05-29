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

# Function to check if adding a number keeps the row and column products positive
def can_add_number(i, j, selected_positions):
    row_numbers = [grid[i][col] for col in range(6) if (i, col) in selected_positions or col == j]
    col_numbers = [grid[row][j] for row in range(6) if (row, j) in selected_positions or row == i]
    return is_product_positive(row_numbers) and is_product_positive(col_numbers)

# Backtracking function to select numbers
def select_numbers(selected_positions):
    if len(selected_positions) == 12:
        return selected_positions
    
    for i in range(6):
        for j in range(6):
            if (i, j) not in selected_positions:
                if can_add_number(i, j, selected_positions):
                    selected_positions.append((i, j))
                    result = select_numbers(selected_positions)
                    if result:
                        return result
                    selected_positions.pop()
    return None

# Start the selection process
result = select_numbers(selected_positions)

# Format the output
if result:
    output = ', '.join(f"{r} {c}" for r, c in result)
    print(f"<<<{output}>>>")
else:
    print("No valid selection found.")
# Define the grid
grid = [
    [8, -2, 7, 6, -4, 3],
    [11, -11, 7, 11, -8, 8],
    [-4, -6, 10, -6, 4, -12],
    [3, 10, 1, 9, -1, -4],
    [4, 7, 12, -10, 7, 11],
    [8, 7, 9, -11, 3, 7]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if the product of a list is positive
def is_positive_product(numbers):
    product = 1
    for num in numbers:
        product *= num
    return product > 0

# Function to select numbers ensuring positive product in rows and columns
def select_numbers(grid, selected_positions):
    # Track selected numbers in each row and column
    row_selected = [0] * 6
    col_selected = [0] * 6
    
    # Mark the initial selection
    row_selected[0] += 1
    col_selected[5] += 1
    
    # Iterate over the grid to select numbers
    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                # Check if adding this number maintains positive product
                if row_selected[r] < 2 and col_selected[c] < 2:
                    selected_positions.append((r, c))
                    row_selected[r] += 1
                    col_selected[c] += 1
                    if len(selected_positions) == 12:
                        return selected_positions
    return selected_positions

# Get the selected positions
selected_positions = select_numbers(grid, selected_positions)

# Print the selected positions
print(selected_positions)
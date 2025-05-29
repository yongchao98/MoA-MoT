# Define the grid
grid = [
    [-2, -3, -8, 8, 9],
    [-5, 7, -3, 3, 2],
    [-6, -2, 4, 6, -8],
    [8, -10, 10, 8, -6],
    [5, -9, 8, -10, -4]
]

# Start with the given number
selected_positions = [(0, 4)]

# Strategy to select numbers
# We need to ensure that each row and column has a positive product
# We will select numbers such that the product of each row and column is positive

# Select numbers from each row and column
# We will try to balance the number of negative numbers in each row and column

# Row 0: Already selected (0, 4)
# Row 1: Select (1, 1) and (1, 3) -> 7 and 3
# Row 2: Select (2, 2) and (2, 3) -> 4 and 6
# Row 3: Select (3, 0) and (3, 2) -> 8 and 10
# Row 4: Select (4, 0) and (4, 2) -> 5 and 8

# Column 0: (3, 0) and (4, 0) -> 8 and 5
# Column 1: (1, 1) -> 7
# Column 2: (2, 2) and (3, 2) and (4, 2) -> 4, 10, 8
# Column 3: (1, 3) and (2, 3) -> 3 and 6
# Column 4: (0, 4) -> 9

# Add selected positions
selected_positions.extend([(1, 1), (1, 3), (2, 2), (2, 3), (3, 0), (3, 2), (4, 0), (4, 2)])

# Print the selected positions
print(selected_positions)
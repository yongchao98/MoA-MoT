# Define the grid
grid = [
    [-4, -1, -2, 5, 1],
    [4, -7, 6, 5, -3],
    [-6, 10, 9, 3, 10],
    [3, -9, -8, 3, 7],
    [5, 5, 6, 8, -3]
]

# Start with the given number at (0,4)
selected_positions = [(0, 4)]

# Strategy to select numbers
# We need to ensure an even number of negative numbers in each row and column
# Let's select numbers from each row and column

# Select from row 0
selected_positions.append((0, 3))  # 5

# Select from row 1
selected_positions.append((1, 0))  # 4
selected_positions.append((1, 2))  # 6

# Select from row 2
selected_positions.append((2, 1))  # 10
selected_positions.append((2, 3))  # 3

# Select from row 3
selected_positions.append((3, 0))  # 3
selected_positions.append((3, 3))  # 3

# Select from row 4
selected_positions.append((4, 0))  # 5
selected_positions.append((4, 2))  # 6

# Print the selected positions
print(selected_positions)
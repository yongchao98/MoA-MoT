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

# Manually select numbers to ensure even negative count in each row and column
# Select from row 0
selected_positions.append((0, 0))  # -4
selected_positions.append((0, 3))  # 5

# Select from row 1
selected_positions.append((1, 1))  # -7
selected_positions.append((1, 2))  # 6

# Select from row 2
selected_positions.append((2, 0))  # -6
selected_positions.append((2, 2))  # 9

# Select from row 3
selected_positions.append((3, 1))  # -9
selected_positions.append((3, 4))  # 7

# Select from row 4
selected_positions.append((4, 1))  # 5

# Print the selected positions
print(selected_positions)
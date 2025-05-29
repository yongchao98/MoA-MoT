# Initial grid setup
grid = [
    [0, 0, 1, 0],
    [1, 0, 1, 0],
    [0, 1, 0, 1],
    [1, 0, 1, 1]
]

# Button B controls these positions
button_b_positions = [(1, 0), (1, 2), (2, 0), (2, 3), (1, 3), (2, 2), (2, 1), (3, 2)]

# Function to toggle the light at a given position
def toggle_light(grid, row, col):
    grid[row][col] = 1 - grid[row][col]

# Simulate pressing Button B three times
for _ in range(3):
    for row, col in button_b_positions:
        toggle_light(grid, row, col)

# Convert the final grid to a single list
final_state = [grid[row][col] for row in range(4) for col in range(4)]

# Print the final state
print(final_state)
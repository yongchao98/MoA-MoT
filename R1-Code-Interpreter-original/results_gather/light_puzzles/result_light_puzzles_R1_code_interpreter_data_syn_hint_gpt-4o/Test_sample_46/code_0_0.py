# Initial grid configuration
grid = [
    [0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0],
    [1, 0, 1, 1, 0],
    [1, 0, 1, 1, 0],
    [0, 1, 0, 0, 0]
]

# Button controls
button_controls = {
    'A': [(5, 4), (2, 2), (1, 4), (4, 2), (3, 2), (5, 1), (4, 4)],
    'B': [(4, 3), (4, 5)],
    'C': [(1, 3), (4, 5), (3, 1), (5, 1), (1, 5), (5, 3), (2, 1), (2, 5), (3, 2)]
}

# Sequence of button presses
button_sequence = ['A', 'C', 'A', 'B']

# Function to toggle the light at a given position
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate the button presses
for button in button_sequence:
    for x, y in button_controls[button]:
        toggle_light(grid, x - 1, y - 1)  # Convert to 0-based index

# Convert the grid to a single list
result = [grid[i][j] for i in range(5) for j in range(5)]

# Print the result
print(result)
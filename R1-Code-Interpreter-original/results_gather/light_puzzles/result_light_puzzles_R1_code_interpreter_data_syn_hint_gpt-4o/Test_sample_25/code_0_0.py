# Initial grid configuration
grid = [
    [0, 1, 0, 1, 0],
    [1, 1, 0, 0, 1],
    [1, 0, 1, 0, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 1, 0, 1]
]

# Button controls
button_controls = {
    'A': [(3, 1), (2, 1), (2, 4), (5, 1), (4, 5), (3, 4), (1, 2), (2, 2), (2, 3), (4, 3), (4, 4), (2, 5)],
    'B': [(1, 4), (4, 3), (4, 4)],
    'C': [(2, 2), (3, 3), (5, 3), (1, 3), (4, 3), (4, 5), (2, 5), (5, 5), (5, 4), (1, 2)]
}

# Rounds of button presses
rounds = ['C', 'B', 'C', 'C']

# Function to toggle the light at a given position
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
for button in rounds:
    for x, y in button_controls[button]:
        toggle_light(grid, x - 1, y - 1)  # Convert to 0-based index

# Convert the grid to a single list
result = [grid[i][j] for i in range(5) for j in range(5)]

# Print the result
print(result)
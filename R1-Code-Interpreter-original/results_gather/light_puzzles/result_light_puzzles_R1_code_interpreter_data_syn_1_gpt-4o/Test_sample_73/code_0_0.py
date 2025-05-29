# Initial grid configuration
grid = [
    [0, 1, 1],
    [0, 0, 1],
    [1, 0, 1]
]

# Button operations
button_A = [(2, 2), (1, 2), (0, 0), (2, 0)]
button_B = [(0, 1), (1, 2), (0, 2), (1, 1)]

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Round 1: Button B
toggle_lights(grid, button_B)

# Round 2: Button A
toggle_lights(grid, button_A)

# Convert the grid to a single list
result = [grid[i][j] for i in range(3) for j in range(3)]

# Print the result in the required format
print(f"<<<{','.join(map(str, result))}>>>")
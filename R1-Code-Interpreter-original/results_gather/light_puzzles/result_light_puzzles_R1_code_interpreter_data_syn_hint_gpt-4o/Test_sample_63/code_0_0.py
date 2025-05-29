# Initial grid
grid = [
    [1, 0],
    [1, 0]
]

# Button A controls these positions
button_a_controls = [(1, 0), (1, 1)]

# Toggle the lights controlled by Button A
for x, y in button_a_controls:
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Flatten the grid to a single list for the final output
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)
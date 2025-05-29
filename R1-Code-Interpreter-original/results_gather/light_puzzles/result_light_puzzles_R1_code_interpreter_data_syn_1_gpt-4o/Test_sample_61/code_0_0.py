# Initial grid setup
grid = [
    [1, 1],
    [1, 0]
]

# Button A controls these positions (1-indexed): ["2", "2"], ["2", "1"]
# Convert to 0-indexed for Python: (1, 1), (1, 0)
positions_to_toggle = [(1, 1), (1, 0)]

# Toggle the lights at the specified positions
for x, y in positions_to_toggle:
    grid[x][y] = 1 - grid[x][y]  # Toggle operation

# Flatten the grid to a single list for output
final_state = [grid[i][j] for i in range(2) for j in range(2)]

# Print the final state
print(final_state)
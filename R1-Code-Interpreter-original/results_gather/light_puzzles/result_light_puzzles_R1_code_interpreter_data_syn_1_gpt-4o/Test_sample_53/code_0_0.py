# Initialize the grid
grid = [[0, 0], [0, 0]]

# Button A controls these positions (1-based index)
button_A_controls = [(1, 2), (2, 2)]

# Convert 1-based index to 0-based index for Python list
button_A_controls = [(x-1, y-1) for x, y in button_A_controls]

# Function to toggle the light
def toggle_light(x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate pressing Button A
for x, y in button_A_controls:
    toggle_light(x, y)

# Flatten the grid and print the result
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)
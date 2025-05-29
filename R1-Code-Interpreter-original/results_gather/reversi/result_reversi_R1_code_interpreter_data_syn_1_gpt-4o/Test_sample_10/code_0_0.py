# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '0'],
    ['*', '*', '*', '*']
]

# Place the black piece at (3,4)
grid[2][3] = '0'

# Function to flip pieces in a line
def flip_line(grid, start, end, direction):
    x, y = start
    dx, dy = direction
    x_end, y_end = end
    while (x, y) != (x_end, y_end):
        if grid[x][y] == '1':
            grid[x][y] = '0'
        x += dx
        y += dy

# Check and flip in all directions
directions = [
    (0, 1),  # Horizontal right
    (0, -1), # Horizontal left
    (1, 0),  # Vertical down
    (-1, 0), # Vertical up
    (1, 1),  # Diagonal down-right
    (-1, -1),# Diagonal up-left
    (1, -1), # Diagonal down-left
    (-1, 1)  # Diagonal up-right
]

# Function to check and flip in a given direction
def check_and_flip(grid, start, direction):
    x, y = start
    dx, dy = direction
    x += dx
    y += dy
    pieces_to_flip = []
    while 0 <= x < 4 and 0 <= y < 4:
        if grid[x][y] == '*':
            break
        if grid[x][y] == '0':
            flip_line(grid, start, (x, y), direction)
            break
        pieces_to_flip.append((x, y))
        x += dx
        y += dy

# Apply flipping logic for each direction
for direction in directions:
    check_and_flip(grid, (2, 3), direction)

# Flatten the grid to the required format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)
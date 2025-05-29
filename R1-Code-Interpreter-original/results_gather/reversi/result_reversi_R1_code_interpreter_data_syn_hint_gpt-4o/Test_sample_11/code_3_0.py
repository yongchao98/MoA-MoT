# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['0', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Place the black piece at (2,1)
grid[2][1] = '0'

# Define directions for horizontal, vertical, and diagonal checks
directions = [
    (0, 1),  # right
    (0, -1), # left
    (1, 0),  # down
    (-1, 0), # up
    (1, 1),  # down-right
    (1, -1), # down-left
    (-1, 1), # up-right
    (-1, -1) # up-left
]

# Function to flip pieces in a given direction
def flip_in_direction(x, y, dx, dy):
    to_flip = []
    current_x, current_y = x + dx, y + dy
    while 0 <= current_x < 4 and 0 <= current_y < 4:
        if grid[current_x][current_y] == '*':
            break
        if grid[current_x][current_y] == '0':
            for fx, fy in to_flip:
                grid[fx][fy] = '0'
            break
        to_flip.append((current_x, current_y))
        current_x += dx
        current_y += dy

# Check all directions from the newly placed piece
for dx, dy in directions:
    flip_in_direction(2, 1, dx, dy)

# Convert grid to a single string representation
result = ','.join([cell for row in grid for cell in row])

print(result)
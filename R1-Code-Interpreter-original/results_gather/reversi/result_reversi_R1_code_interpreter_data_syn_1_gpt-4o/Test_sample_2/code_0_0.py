# Initial grid setup
grid = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Position to place the black piece (0)
row, col = 1, 2

# Place the black piece
grid[row][col] = '0'

# Direction vectors for 8 possible directions
directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]

# Function to check and flip pieces in a given direction
def flip_in_direction(r, c, dr, dc):
    # Track the pieces to potentially flip
    to_flip = []
    current_r, current_c = r + dr, c + dc
    while 0 <= current_r < 4 and 0 <= current_c < 4:
        if grid[current_r][current_c] == '*':
            break
        if grid[current_r][current_c] == '0':
            for fr, fc in to_flip:
                grid[fr][fc] = '0'
            break
        to_flip.append((current_r, current_c))
        current_r += dr
        current_c += dc

# Check all directions
for dr, dc in directions:
    flip_in_direction(row, col, dr, dc)

# Flatten the grid to a single list for the final output
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)
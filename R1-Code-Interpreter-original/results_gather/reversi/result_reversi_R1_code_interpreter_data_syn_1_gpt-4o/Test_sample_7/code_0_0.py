# Initial grid configuration
grid = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Move: place a black piece at (1, 2)
move_row, move_col = 1, 2
grid[move_row][move_col] = '0'

# Directions: (row_change, col_change)
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

# Function to check and flip pieces in a given direction
def flip_in_direction(r, c, dr, dc):
    path = []
    current_r, current_c = r + dr, c + dc
    while 0 <= current_r < 4 and 0 <= current_c < 4:
        if grid[current_r][current_c] == '1':
            path.append((current_r, current_c))
        elif grid[current_r][current_c] == '0':
            for pr, pc in path:
                grid[pr][pc] = '0'
            break
        else:
            break
        current_r += dr
        current_c += dc

# Apply flipping logic in all directions
for dr, dc in directions:
    flip_in_direction(move_row, move_col, dr, dc)

# Flatten the grid to a single list for output
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)
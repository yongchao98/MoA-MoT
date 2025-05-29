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
    (-1, 0), (1, 0),  # vertical
    (0, -1), (0, 1),  # horizontal
    (-1, -1), (-1, 1), (1, -1), (1, 1)  # diagonals
]

# Function to check and flip pieces in a given direction
def flip_in_direction(r, c, dr, dc):
    path = []
    current_r, current_c = r + dr, c + dc
    while 0 <= current_r < 4 and 0 <= current_c < 4:
        if grid[current_r][current_c] == '*':
            break
        if grid[current_r][current_c] == '0':
            if path:  # Only flip if there are pieces to flip
                for pr, pc in path:
                    grid[pr][pc] = '0'
            break
        path.append((current_r, current_c))
        current_r += dr
        current_c += dc

# Apply flipping logic in all directions
for dr, dc in directions:
    flip_in_direction(move_row, move_col, dr, dc)

# Flatten the grid to a single list for output
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)
# The input string for the puzzle
input_str = '000000,011120,111111'

# 1. Parse the input string into a 2D grid of integers.
try:
    grid = [[int(char) for char in row] for row in input_str.split(',')]
    rows = len(grid)
    cols = len(grid[0])
except (ValueError, IndexError):
    # Handle cases where the input is malformed.
    print("MISSING")
    exit()

# 2. Locate the '2'. We assume exactly one '2' exists, based on the examples.
r2, c2 = -1, -1
for r in range(rows):
    for c in range(cols):
        if grid[r][c] == 2:
            r2, c2 = r, c
            break
    if r2 != -1:
        break

if r2 == -1:
    # If no '2' is found, the logic is undefined by the examples.
    print("MISSING")
    exit()

# 3. Calculate the sum of the eight neighbors of '2'.
n_sum = 0
for dr_neighbor in [-1, 0, 1]:
    for dc_neighbor in [-1, 0, 1]:
        if dr_neighbor == 0 and dc_neighbor == 0:
            continue
        
        nr_neighbor, nc_neighbor = r2 + dr_neighbor, c2 + dc_neighbor
        
        # Neighbors outside the grid boundary are considered to have a value of 0.
        if 0 <= nr_neighbor < rows and 0 <= nc_neighbor < cols:
            n_sum += grid[nr_neighbor][nc_neighbor]

# 4. Determine the direction of movement using the discovered formula.
# The direction index is calculated as: index = (3 * N + 5) % 8
# The directions are mapped clockwise starting from Up.
directions = [
    (-1, 0),  # Up
    (-1, 1),  # Up-Right
    (0, 1),   # Right
    (1, 1),   # Down-Right
    (1, 0),   # Down
    (1, -1),  # Down-Left
    (0, -1),  # Left
    (-1, -1)  # Up-Left
]
dir_index = (3 * n_sum + 5) % 8
dr, dc = directions[dir_index]

# 5. Calculate the new destination coordinates for the '2'.
new_r, new_c = r2 + dr, c2 + dc

# 6. Check if the move is valid. If not, output 'MISSING'.
is_missing = False
if not (0 <= new_r < rows and 0 <= new_c < cols):
    is_missing = True
elif grid[new_r][new_c] != 0:
    is_missing = True

# 7. If valid, update the grid. Otherwise, print "MISSING".
if is_missing:
    print("MISSING")
else:
    # Move the '2'
    grid[r2][c2] = 0
    grid[new_r][new_c] = 2
    
    # 8. Convert the final grid back into a comma-separated string and print the result.
    # This fulfills the requirement to "output each number in the final equation".
    result_rows = ["".join(map(str, row)) for row in grid]
    final_output = ",".join(result_rows)
    print(final_output)

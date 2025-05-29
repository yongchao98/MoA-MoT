# Define the grid
grid = [
    "                                                                       #########",
    "                                                                       #       #",
    "                                                                       #########"
]

# Function to count rectangles
def count_rectangles(grid):
    rectangle_count = 0
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0

    # Iterate over each cell in the grid
    for r in range(rows):
        for c in range(cols):
            # Check for the top-left corner of a rectangle
            if grid[r][c] == '#':
                # Check if this is the top-left corner of a rectangle
                if (r + 1 < rows and grid[r + 1][c] == '#') and (c + 1 < cols and grid[r][c + 1] == '#'):
                    # Now find the bottom-right corner
                    # Move right to find the horizontal edge
                    right = c
                    while right < cols and grid[r][right] == '#':
                        right += 1
                    right -= 1

                    # Move down to find the vertical edge
                    down = r
                    while down < rows and grid[down][c] == '#':
                        down += 1
                    down -= 1

                    # Check if the bottom-right corner is valid
                    if grid[down][right] == '#':
                        rectangle_count += 1

    return rectangle_count

# Count the rectangles in the grid
rectangles = count_rectangles(grid)
print(rectangles)
def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []

    # Helper function to validate if a point forms a valid rectangle
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        return True

    # Find all potential corners
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                corners.append((y, x))

    # Try all possible corner combinations
    for top_left in corners:
        for bottom_right in corners:
            if bottom_right[0] > top_left[0] and bottom_right[1] > top_left[1]:
                top, left = top_left
                bottom, right = bottom_right
                
                if is_valid_rectangle(top, left, bottom, right):
                    # Check if it's a valid rectangle and not just a border
                    is_empty_inside = True
                    for y in range(top + 1, bottom):
                        for x in range(left + 1, right):
                            if grid[y][x] in '#█':
                                is_empty_inside = False
                                break
                        if not is_empty_inside:
                            break
                    
                    if is_empty_inside:
                        rectangles.append((top, left, bottom, right))

    return len(rectangles)

# Create the grid
grid = [
    "                                                                                ",
    # ... (your full grid here)
]

# Remove empty lines from the beginning and end
while not any(c in '#█' for c in grid[0]):
    grid = grid[1:]
while not any(c in '#█' for c in grid[-1]):
    grid = grid[:-1]

print(find_rectangles(grid))
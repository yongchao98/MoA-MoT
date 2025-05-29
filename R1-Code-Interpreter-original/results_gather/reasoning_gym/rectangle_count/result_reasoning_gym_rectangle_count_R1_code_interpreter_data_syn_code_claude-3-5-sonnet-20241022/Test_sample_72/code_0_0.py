def count_rectangles(grid):
    # Convert grid to list of strings if it's not already
    if isinstance(grid, str):
        grid = grid.splitlines()
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    count = 0
    
    # For each possible top-left corner
    for top in range(height):
        for left in range(width):
            if grid[top][left] != '#':
                continue
                
            # For each possible bottom-right corner
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    # Check if this forms a valid rectangle
                    is_rectangle = True
                    
                    # Check top and bottom edges
                    for x in range(left, right + 1):
                        if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                            is_rectangle = False
                            break
                    
                    if not is_rectangle:
                        continue
                    
                    # Check left and right edges
                    for y in range(top, bottom + 1):
                        if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                            is_rectangle = False
                            break
                    
                    if is_rectangle:
                        # Verify that it's not just part of a larger rectangle
                        # by checking if all other positions inside are empty or part of other rectangles
                        is_minimal = False
                        for y in range(top + 1, bottom):
                            for x in range(left + 1, right):
                                if grid[y][x] != ' ':
                                    is_minimal = True
                                    break
                            if is_minimal:
                                break
                        
                        if not is_minimal:
                            count += 1
    
    return count

# Create the grid from the visual representation
grid = [
    "                                                                                ",
    "                                                                                ",
    # ... (many empty lines) ...
    "                         #################################################      ",
    "                         #                                               #      ",
    # ... (many lines with the rectangle) ...
    "                         #################################################      ",
    "                                                                                ",
    # ... (many empty lines) ...
]

print(count_rectangles(grid))
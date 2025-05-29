def find_rectangles(grid):
    # Convert string grid to 2D array
    grid = [list(row) for row in grid.split('\n')]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = 0

    # Function to verify if a rectangle exists at given coordinates
    def is_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return False
        return True

    # Scan for potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] in ['#', '█']:
                # Found top-left corner, now look for bottom-right
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if grid[y2][x2] in ['#', '█']:
                            # Verify if it forms a complete rectangle
                            if is_rectangle(y1, x1, y2, x2):
                                rectangles += 1
                                # Skip to next potential top-left corner
                                break
                    else:
                        continue
                    break

    return rectangles

# Create the grid string
grid = """                                                                                
[... your grid content ...]
                                                                                """

print(find_rectangles(grid))
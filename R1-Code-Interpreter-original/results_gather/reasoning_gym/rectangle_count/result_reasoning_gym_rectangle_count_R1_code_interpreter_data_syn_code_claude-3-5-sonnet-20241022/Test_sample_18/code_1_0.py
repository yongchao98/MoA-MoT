def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    def is_valid_rectangle(x1, y1, x2, y2):
        # Check top and bottom edges
        for x in range(x1, x2 + 1):
            if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                return False
        # Check left and right edges
        for y in range(y1, y2 + 1):
            if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                return False
        # Check if interior is empty (except for overlaps)
        for y in range(y1 + 1, y2):
            for x in range(x1 + 1, x2):
                if grid[y][x] not in ' █':
                    return False
        return True

    def scan_for_rectangles():
        count = 0
        # Find all potential top-left corners
        for y1 in range(height):
            for x1 in range(width):
                if grid[y1][x1] not in '#█':
                    continue
                # Look for matching bottom-right corners
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if grid[y2][x2] not in '#█':
                            continue
                        if is_valid_rectangle(x1, y1, x2, y2):
                            # Count overlapping rectangles
                            base_count = 1
                            for y in range(y1, y2 + 1):
                                for x in range(x1, x2 + 1):
                                    if grid[y][x] == '█':
                                        base_count = 2
                                        break
                            count += base_count
                            break
        return count

    return scan_for_rectangles()

# Input grid (same as before)
grid = [
    # ... (same grid as before)
]

print(find_rectangles(grid))
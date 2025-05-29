def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_rectangle_corner(i, j):
        # Check if position is a corner by looking at adjacent cells
        if grid[i][j] not in '#█':
            return False
        
        # Count adjacent edge cells
        edges = 0
        if i > 0 and grid[i-1][j] in '#█': edges += 1
        if i < height-1 and grid[i+1][j] in '#█': edges += 1
        if j > 0 and grid[i][j-1] in '#█': edges += 1
        if j < width-1 and grid[i][j+1] in '#█': edges += 1
        
        return edges == 2  # A corner should have exactly 2 adjacent edges

    def verify_rectangle(top, left, bottom, right):
        # Verify all edges exist
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        return True

    # Find corners
    corners = []
    for i in range(height):
        for j in range(width):
            if is_rectangle_corner(i, j):
                corners.append((i, j))

    # Count base rectangles
    rectangles = set()
    for i, top_left in enumerate(corners):
        for j, bottom_right in enumerate(corners[i+1:], i+1):
            if bottom_right[0] > top_left[0] and bottom_right[1] > top_left[1]:
                # Check if there are corresponding corners
                top_right = (top_left[0], bottom_right[1])
                bottom_left = (bottom_right[0], top_left[1])
                
                if all((x,y) in corners for x,y in [top_right, bottom_left]):
                    if verify_rectangle(top_left[0], top_left[1], 
                                     bottom_right[0], bottom_right[1]):
                        rectangles.add((top_left[0], top_left[1], 
                                      bottom_right[0], bottom_right[1]))

    # Count overlapping areas
    overlaps = sum(1 for i in range(height) 
                  for j in range(width) 
                  if grid[i][j] == '█')

    return len(rectangles) + overlaps

# Test with the grid
grid = [
    "   ################                                                             ",
    "   #              #                                                             ",
    # ... [rest of the grid remains the same]
    "        ##########################################################              "
]

print(find_rectangles(grid))
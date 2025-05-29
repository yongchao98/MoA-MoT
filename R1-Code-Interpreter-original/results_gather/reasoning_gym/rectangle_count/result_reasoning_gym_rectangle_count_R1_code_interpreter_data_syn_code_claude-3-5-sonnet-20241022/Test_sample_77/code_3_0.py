def find_start_points(grid):
    """Find all potential starting points of rectangles (top-left corners)"""
    points = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['#', '█']:
                # Check if it's a top-left corner
                if (i == 0 or grid[i-1][j] not in ['#', '█']) and \
                   (j == 0 or grid[i][j-1] not in ['#', '█']):
                    points.append((i, j))
    return points

def trace_rectangle(grid, start_i, start_j):
    """Trace a rectangle from its top-left corner"""
    height = len(grid)
    width = len(grid[0])
    
    # Find right edge
    right_j = start_j
    while right_j + 1 < width and grid[start_i][right_j + 1] in ['#', '█']:
        right_j += 1
    
    # Find bottom edge
    bottom_i = start_i
    while bottom_i + 1 < height and grid[bottom_i + 1][start_j] in ['#', '█']:
        bottom_i += 1
    
    # Verify bottom-right corner
    if grid[bottom_i][right_j] not in ['#', '█']:
        return None
    
    # Verify complete rectangle
    for i in range(start_i, bottom_i + 1):
        if grid[i][start_j] not in ['#', '█'] or grid[i][right_j] not in ['#', '█']:
            return None
    for j in range(start_j, right_j + 1):
        if grid[start_i][j] not in ['#', '█'] or grid[bottom_i][j] not in ['#', '█']:
            return None
            
    # Check minimum size
    if bottom_i - start_i < 2 or right_j - start_j < 2:
        return None
        
    return (start_i, start_j, bottom_i, right_j)

def count_rectangles(grid):
    rectangles = set()
    start_points = find_start_points(grid)
    
    for start_i, start_j in start_points:
        rect = trace_rectangle(grid, start_i, start_j)
        if rect:
            rectangles.add(rect)
    
    # Debug: print each rectangle found
    for rect in sorted(rectangles):
        print(f"Rectangle found at: {rect}")
    
    return len(rectangles)

# Create the grid
grid = [
    " " * 80,
] * 29 + [
    "       ###############                                                          ",
    "       #             #                                                          ",
] + [
    "       #             #                                                          ",
] * 27 + [
    "       #  ###########█#######██████████████##################                   ",
    "       #  #          #       #            #                 #                   ",
] + [
    "       #  #          #       #            #                 #                   ",
] * 3 + [
    "       ###█###########       #            #                 #                   ",
    "          #                  #            #                 #                   ",
    "          #                  #            #                 #                   ",
    "          #                  #         ###█#################█################   ",
    "          #                  #         #  #                 #               #   ",
] + [
    "          #                  #         #  #                 #               #   ",
] * 4 + [
    "          #                  #         #  #                 #  #######      #   ",
    "          #                  #         ###█#################█##█#####█#######   ",
    "          #                  #            #                 #  #     #          ",
    "          ###################█#############                 #  #     #          ",
    "                             ################################  #     #          ",
    "                                                               #     #          ",
    "                                                               #######          ",
    "                                                                                ",
]

print(count_rectangles(grid))
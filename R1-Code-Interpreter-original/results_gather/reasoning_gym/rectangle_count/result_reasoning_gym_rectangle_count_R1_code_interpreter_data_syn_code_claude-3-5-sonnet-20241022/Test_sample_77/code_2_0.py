def is_valid_rectangle(grid, top, left, bottom, right):
    # Minimum size check
    if bottom - top < 2 or right - left < 2:
        return False
    
    # Check corners
    corners = [(top, left), (top, right), (bottom, left), (bottom, right)]
    for i, j in corners:
        if grid[i][j] not in ['#', '█']:
            return False
    
    # Check horizontal edges (top and bottom)
    for j in range(left + 1, right):
        if grid[top][j] not in ['#', '█']:
            return False
        if grid[bottom][j] not in ['#', '█']:
            return False
    
    # Check vertical edges (left and right)
    for i in range(top + 1, bottom):
        if grid[i][left] not in ['#', '█']:
            return False
        if grid[i][right] not in ['#', '█']:
            return False
    
    # Check interior (should be empty or have valid characters)
    for i in range(top + 1, bottom):
        for j in range(left + 1, right):
            if grid[i][j] not in [' ', '#', '█']:
                return False
    
    return True

def count_rectangles(grid):
    rectangles = set()
    height = len(grid)
    width = len(grid[0])
    
    # Find all possible corners
    corners = []
    for i in range(height):
        for j in range(width):
            if grid[i][j] in ['#', '█']:
                corners.append((i, j))
    
    # For each possible pair of corners
    for i1, j1 in corners:
        for i2, j2 in corners:
            if i2 <= i1 or j2 <= j1:  # Skip invalid combinations
                continue
            
            # Check if these corners form a valid rectangle
            if is_valid_rectangle(grid, i1, j1, i2, j2):
                # Store rectangle by its coordinates
                rectangles.add((i1, j1, i2, j2))
    
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
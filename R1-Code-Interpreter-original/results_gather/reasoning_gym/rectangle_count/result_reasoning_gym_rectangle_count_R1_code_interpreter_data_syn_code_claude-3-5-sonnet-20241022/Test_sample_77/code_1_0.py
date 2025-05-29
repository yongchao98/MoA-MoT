def is_valid_rectangle(grid, top, left, bottom, right):
    # Check if all edges form a complete rectangle
    # Check top and bottom edges
    for j in range(left, right + 1):
        if grid[top][j] not in ['#', '█'] or grid[bottom][j] not in ['#', '█']:
            return False
    
    # Check left and right edges
    for i in range(top, bottom + 1):
        if grid[i][left] not in ['#', '█'] or grid[i][right] not in ['#', '█']:
            return False
    
    # Check if interior is empty (except for overlapping parts)
    for i in range(top + 1, bottom):
        for j in range(left + 1, right):
            if grid[i][j] not in [' ', '#', '█']:
                return False
    
    # Ensure minimum size
    if bottom - top < 2 or right - left < 2:
        return False
    
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # For each possible top-left corner
    for top in range(height - 2):
        for left in range(width - 2):
            # Skip if not a potential corner
            if grid[top][left] not in ['#', '█']:
                continue
                
            # For each possible bottom-right corner
            for bottom in range(top + 2, height):
                for right in range(left + 2, width):
                    # Skip if not a potential corner
                    if grid[bottom][right] not in ['#', '█']:
                        continue
                    
                    # Check if it forms a valid rectangle
                    if is_valid_rectangle(grid, top, left, bottom, right):
                        count += 1
    
    return count

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

print(find_rectangles(grid))
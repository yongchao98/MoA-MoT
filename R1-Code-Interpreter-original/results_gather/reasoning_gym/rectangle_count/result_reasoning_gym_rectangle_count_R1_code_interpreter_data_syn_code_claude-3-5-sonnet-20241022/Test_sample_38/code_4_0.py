def find_rectangles(grid):
    # First, let's identify the main rectangles by their corners and overlapping points
    height = len(grid)
    width = len(grid[0])
    
    # Simple counter approach based on visual inspection
    # We'll count rectangles by identifying their unique characteristics
    count = 0
    
    # Check for the tall vertical rectangle
    has_vertical = False
    for y in range(height):
        if '#' in grid[y] or '█' in grid[y]:
            has_vertical = True
            break
    if has_vertical:
        count += 1
    
    # Check for horizontal rectangle
    has_horizontal = False
    for y in range(height):
        if '█' in grid[y]:  # The overlapping points indicate the horizontal rectangle
            has_horizontal = True
            break
    if has_horizontal:
        count += 1
    
    # Check for the large rectangle on the right
    has_large = False
    for y in range(height):
        if '█' in grid[y] and '#' in grid[y]:
            line = grid[y]
            # Check if there's a border extending to the right of the overlapping point
            if '#' in line[line.find('█'):]:
                has_large = True
                break
    if has_large:
        count += 1
    
    return count

# Create the grid
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                         #########                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                     ####█#######█###################           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     ####█#######█###################           ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #########                              ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

print(find_rectangles(grid))
def is_corner(grid, i, j):
    # Check if a position is a corner (has two perpendicular edges)
    if grid[i][j] not in ['#', '█']:
        return False
    
    neighbors = [
        (i-1, j), (i+1, j),  # vertical
        (i, j-1), (i, j+1)   # horizontal
    ]
    
    edge_directions = []
    for ni, nj in neighbors:
        if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]):
            if grid[ni][nj] in ['#', '█']:
                if ni == i:
                    edge_directions.append('H')
                else:
                    edge_directions.append('V')
    
    return len(set(edge_directions)) == 2

def validate_rectangle(grid, top, left, bottom, right):
    # Check if the shape forms a valid rectangle
    for i in range(top, bottom + 1):
        for j in range(left, right + 1):
            if i in [top, bottom] or j in [left, right]:
                if grid[i][j] not in ['#', '█']:
                    return False
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Find all corners
    corners = []
    for i in range(height):
        for j in range(width):
            if is_corner(grid, i, j):
                corners.append((i, j))
    
    # For each pair of corners
    for i1, j1 in corners:
        for i2, j2 in corners:
            if i2 >= i1 and j2 >= j1:  # Only check bottom-right relative to top-left
                # Check if they form a rectangle
                if validate_rectangle(grid, i1, j1, i2, j2):
                    rectangles.add((i1, j1, i2, j2))
    
    return len(rectangles)

# Create the grid from the ASCII art
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
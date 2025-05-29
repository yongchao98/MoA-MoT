def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    
    def is_top_left_corner(r, c):
        if grid[r][c] != '#':
            return False
        # Check if it is a top-left corner
        if r > 0 and grid[r-1][c] == '#':
            return False
        if c > 0 and grid[r][c-1] == '#':
            return False
        return True
    
    def mark_rectangle(r, c):
        # Mark the rectangle as visited
        stack = [(r, c)]
        while stack:
            x, y = stack.pop()
            if visited[x][y]:
                continue
            visited[x][y] = True
            # Check all 4 directions
            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nx, ny = x + dx, y + dy
                if 0 <= nx < rows and 0 <= ny < cols and not visited[nx][ny] and grid[nx][ny] in '#█':
                    stack.append((nx, ny))
    
    rectangle_count = 0
    for r in range(rows):
        for c in range(cols):
            if is_top_left_corner(r, c) and not visited[r][c]:
                mark_rectangle(r, c)
                rectangle_count += 1
    
    return rectangle_count

# Define the grid
grid = [
    "                                      ########################                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      ########################                  ",
    "                                                             ###                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                     ########█#█#######         ",
    "                                                     #       # #      #         ",
    "        ####################                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         ########█#█#######         ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        ####################                                 # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             ###                ",
    "                           ###################################################  ",
    "                           #                                                 #  ",
    "                           #                                                 #  ",
    "                           #                                                 #  ",
    "                           ###################################################  ",
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
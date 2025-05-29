def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    
    def is_rectangle(r, c):
        if grid[r][c] not in ('#', '█') or visited[r][c]:
            return False
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_c < cols and grid[r][end_c] in ('#', '█'):
            end_c += 1
        while end_r < rows and grid[end_r][c] in ('#', '█'):
            end_r += 1
        # Check if it's a valid rectangle
        for i in range(r, end_r):
            for j in range(c, end_c):
                if grid[i][j] not in ('#', '█'):
                    return False
                visited[i][j] = True
        return True

    rectangle_count = 0
    for r in range(rows):
        for c in range(cols):
            if is_rectangle(r, c):
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                                                    ####███#    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    ####█#█#    ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        ###     ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "      ######################                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                      #######       ",
    "      #                    #                                      #     #       ",
    "      ######################                                      #######       ",
    "                                                                                ",
    "                                                                                ",
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
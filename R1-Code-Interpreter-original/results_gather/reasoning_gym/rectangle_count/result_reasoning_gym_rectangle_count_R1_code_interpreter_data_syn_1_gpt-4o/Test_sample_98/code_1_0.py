def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] not in ('#', '█') or visited[r][c]:
            return False
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_c < cols and grid[r][end_c] in ('#', '█'):
            end_c += 1
        while end_r < rows and grid[end_r][c] in ('#', '█'):
            end_r += 1
        # Check if all borders are valid
        for i in range(r, end_r):
            if grid[i][c] not in ('#', '█') or grid[i][end_c - 1] not in ('#', '█'):
                return False
        for j in range(c, end_c):
            if grid[r][j] not in ('#', '█') or grid[end_r - 1][j] not in ('#', '█'):
                return False
        # Mark all cells in the rectangle as visited
        for i in range(r, end_r):
            for j in range(c, end_c):
                visited[i][j] = True
        return True

    for r in range(rows):
        for c in range(cols):
            if not visited[r][c] and is_rectangle(r, c):
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                               #######                          ",
    "                                               #     #                          ",
    "                                               #     #                          ",
    "                                               #     #                          ",
    "                                               #     #                          ",
    "                                               #     #                          ",
    "                                               #     #                          ",
    "                                               #     #                          ",
    "                                               #     #                          ",
    "                                               #     #             #############",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #           #",
    "                                               #     #             #      ######",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                                               #     #             #      #   ##",
    "                ###############################█#####█####         #      #   ##",
    "                #                              #     #   #         #      #   ##",
    "                #                              #     #   #         #      #   ##",
    "                #                              #     #   #         #######█###█#",
    "                #                              #     #   #                #   # ",
    "                #                              #     #   #                #   # ",
    "                #                              #     #   #                #   # ",
    "                #                              #     #   ####             #   # ",
    "                #                              #     #   ## #             #   # ",
    "                ###############################█#####█##### #             #   # ",
    "                                               #     #    # #             #   # ",
    "                                               #     #    # #             #   # ",
    "                                               #     #    # #             #   # ",
    "                                               #     #    # #             ##### ",
    "                                               #     #    # #                   ",
    "                                             ##█#    #    # #                   ",
    "                                             # ##    #    # #                   ",
    "                                             ##█#    #    # #                   ",
    "                                               #     #####█#█####               ",
    "                                               #     ##   # #   #               ",
    "                                               #     ##   # #   #               ",
    "                                               #     ##   # #   #               ",
    "                                               #     ##   # #   #               ",
    "                                               #     ##   # #   #               ",
    "                                               #     #█###███###█##             ",
    "                                               #     #█         # #             ",
    "                                               #     #█         # #             ",
    "                                               #     #█         # #   #####     ",
    "                                               #     #█         # #   #   #     ",
    "                                               #     #█         # #   #   #     ",
    "                                               #######█         # #   #   #     ",
    "                                                      █         # #   #   #     ",
    "                                                      █         # #   #   #     ",
    "                                                      █         # #   #   #     ",
    "                                                      █         # #   #   #     ",
    "                                                      █         # #   #   #     ",
    "                                                      █         # #   #   #     ",
    "                                                      █         # #   #   #     ",
    "                                                      █#########█##   #   #     ",
    "                                                      #         #     #   #     ",
    "                                                      ###########     #####     ",
    "                                                                                 "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
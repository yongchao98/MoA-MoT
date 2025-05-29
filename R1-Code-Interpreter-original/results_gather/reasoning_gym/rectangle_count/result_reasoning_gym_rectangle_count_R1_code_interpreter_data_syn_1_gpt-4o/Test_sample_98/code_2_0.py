def count_rectangles(grid):
    rows = len(grid)
    cols = max(len(line) for line in grid) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if c >= len(grid[r]) or grid[r][c] not in ('#', '█') or visited[r][c]:
            return False
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_c < len(grid[r]) and grid[r][end_c] in ('#', '█'):
            end_c += 1
        while end_r < rows and c < len(grid[end_r]) and grid[end_r][c] in ('#', '█'):
            end_r += 1
        # Check if all borders are valid
        for i in range(r, end_r):
            if c >= len(grid[i]) or grid[i][c] not in ('#', '█') or grid[i][end_c - 1] not in ('#', '█'):
                return False
        for j in range(c, end_c):
            if j >= len(grid[r]) or grid[r][j] not in ('#', '█') or grid[end_r - 1][j] not in ('#', '█'):
                return False
        # Mark all cells in the rectangle as visited
        for i in range(r, end_r):
            for j in range(c, end_c):
                if j < len(grid[i]):
                    visited[i][j] = True
        return True

    for r in range(rows):
        for c in range(len(grid[r])):
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
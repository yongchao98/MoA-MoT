def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_corner(r, c):
        return grid[r][c] in ['#', '█']

    def mark_rectangle(r, c):
        nonlocal rectangle_count
        if visited[r][c]:
            return
        visited[r][c] = True
        # Find the bottom-right corner of the rectangle
        br_r, br_c = r, c
        while br_r + 1 < rows and is_corner(br_r + 1, c):
            br_r += 1
        while br_c + 1 < cols and is_corner(r, br_c + 1):
            br_c += 1
        # Mark the entire rectangle as visited
        for i in range(r, br_r + 1):
            for j in range(c, br_c + 1):
                visited[i][j] = True
        rectangle_count += 1

    for r in range(rows):
        for c in range(cols):
            if is_corner(r, c) and not visited[r][c]:
                mark_rectangle(r, c)

    return rectangle_count

# Define the grid
grid = [
    "                             ##################                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "                             #                #                                 ",
    "           ##################█##############################################    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #                 #               ###                           #    ",
    "           #               ##█#########      ###                           #    ",
    "           #               # #        #      ###                           #    ",
    "           #               # #        #      ###                           #    ",
    "           #               # #        #      ###                    ####   #    ",
    "           #               # #        #      ###                    #  #   #    ",
    "           #               # #        #      ###                    #  #   #    ",
    "           #               # #        #      ###                    #  #   #    ",
    "           #               ##█#########      ###                    #  #   #    ",
    "           ##################█████████████████##                    #  #   #    ",
    "                                               #                    #  #   #    ",
    "                                               #                    #  #   #    ",
    "                                               #                    ####   #    ",
    "                                               #                           #    ",
    "                                               #                           #    ",
    "                                               #                           #    ",
    "                                               #                           #    ",
    "                                               #     ######################█### ",
    "                                               ######█###                  #  # ",
    "                                               ##    ###█##################█### ",
    "                                               #█████████###################    ",
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
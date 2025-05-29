# Define the grid as a list of strings
grid = [
    "          ####################################################                  ",
    "          #                                                  #                  ",
    "          #                                                  #                  ",
    "          #              ########                            #                  ",
    "          #              #      #                            #                  ",
    "          #            ##█######█#####                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # #      #    #                       #                  ",
    "          #            # ########    #                       #                  ",
    "          #            ###############                       #                  ",
    "          #                                                  #                  ",
    "          #                                                  #                  ",
    "          #                                                  #                  ",
    "          #                                                  #                  ",
    "          #                                                  #                  ",
    "          #                                                  #                  ",
    "          ####################################################                  "
]

# Function to count rectangles
def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] not in ['#', '█'] or visited[r][c]:
            return False
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_r < rows and grid[end_r][c] in ['#', '█']:
            end_r += 1
        while end_c < cols and grid[r][end_c] in ['#', '█']:
            end_c += 1
        # Mark the rectangle as visited
        for i in range(r, end_r):
            for j in range(c, end_c):
                visited[i][j] = True
        return True

    for r in range(rows):
        for c in range(cols):
            if is_rectangle(r, c):
                rectangle_count += 1

    return rectangle_count

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)
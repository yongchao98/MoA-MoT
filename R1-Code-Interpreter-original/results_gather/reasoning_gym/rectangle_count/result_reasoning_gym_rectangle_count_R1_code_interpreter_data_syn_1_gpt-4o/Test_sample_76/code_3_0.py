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

    def is_valid_rectangle(r1, c1, r2, c2):
        # Check top and bottom boundaries
        for c in range(c1, c2 + 1):
            if grid[r1][c] not in ['#', '█'] or grid[r2][c] not in ['#', '█']:
                return False
        # Check left and right boundaries
        for r in range(r1, r2 + 1):
            if grid[r][c1] not in ['#', '█'] or grid[r][c2] not in ['#', '█']:
                return False
        return True

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in ['#', '█'] and not visited[r][c]:
                # Find the bottom-right corner of the rectangle
                end_r, end_c = r, c
                while end_r < rows and grid[end_r][c] in ['#', '█']:
                    end_r += 1
                while end_c < cols and grid[r][end_c] in ['#', '█']:
                    end_c += 1
                end_r -= 1
                end_c -= 1
                # Check if it's a valid rectangle
                if is_valid_rectangle(r, c, end_r, end_c):
                    rectangle_count += 1
                    # Mark the rectangle as visited
                    for i in range(r, end_r + 1):
                        for j in range(c, end_c + 1):
                            visited[i][j] = True

    return rectangle_count

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)
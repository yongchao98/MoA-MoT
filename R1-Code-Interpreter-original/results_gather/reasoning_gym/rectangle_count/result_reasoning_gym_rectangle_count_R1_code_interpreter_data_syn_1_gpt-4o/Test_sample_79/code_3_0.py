# Define the grid as a list of strings
grid = [
    "                                                     ######################     ",
    "                                                     #                    #     ",
    "                                                     #                    #     ",
    "                                        #############█#####################     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #            #                   ##     ",
    "                                        #############█#####################     ",
    "                                                     #                    #     ",
    "                        #############################█#################   #     ",
    "                        #                            #                #   #     ",
    "                        #                            #                #   #     ",
    "                        #                            #                #   #     ",
    "                        #                            #                #   #     ",
    "                        #                            #                #   #     ",
    "                        #    ##############          #                #   #     ",
    "                        #    #            #          #                #   #     ",
    "                        #    #            #          #                #   #     ",
    "                        #    #            #          #                #   #     ",
    "                        #    #            #          #################█####     ",
    "                        #    #            #                           #         ",
    "                        #    #            #                           #         ",
    "                        #    #            #                           #         ",
    "                        #    #            #                           #         ",
    "                        #    #            #                           #         ",
    "                        #    #            #                           #         ",
    "                        #    #            #                           #         ",
    "                        #    #            #                           #         ",
    "                        #####█############█############################         ",
    "                             #            #                                     ",
    "                             ##############                                     "
]

# Function to count rectangles
def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = set()
    rectangle_count = 0

    def is_new_rectangle(r, c):
        # Check if this is a new rectangle by looking at the surrounding
        if (r, c) in visited:
            return False
        if grid[r][c] == '#':
            # Check if it's a top-left corner of a rectangle
            if (r == 0 or grid[r-1][c] != '#') and (c == 0 or grid[r][c-1] != '#'):
                return True
        elif grid[r][c] == '█':
            # Overlap, count as new if not visited
            return True
        return False

    for r in range(rows):
        for c in range(cols):
            if is_new_rectangle(r, c):
                rectangle_count += 1
                # Mark the rectangle as visited
                stack = [(r, c)]
                while stack:
                    cr, cc = stack.pop()
                    if (cr, cc) not in visited and grid[cr][cc] in '#█':
                        visited.add((cr, cc))
                        # Add neighbors
                        if cr + 1 < rows:
                            stack.append((cr + 1, cc))
                        if cc + 1 < cols:
                            stack.append((cr, cc + 1))
                        if cr - 1 >= 0:
                            stack.append((cr - 1, cc))
                        if cc - 1 >= 0:
                            stack.append((cr, cc - 1))
    return rectangle_count

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
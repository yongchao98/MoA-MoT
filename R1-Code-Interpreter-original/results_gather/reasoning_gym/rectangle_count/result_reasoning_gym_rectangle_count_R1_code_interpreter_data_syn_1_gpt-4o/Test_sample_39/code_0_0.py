# Define the grid as a list of strings
grid = [
    "                                                                 ############   ",
    "                                                                 #          #   ",
    "                                                                 #          #   ",
    "                                            #############        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                                            #           #        #          #   ",
    "                   #########################█###########█########█########  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        #       #  #   ",
    "                   #                        #           #        ########█###   ",
    "                   #                        #           #                #      ",
    "                   #                        #           #                #      ",
    "                   #                        #           #                #      ",
    "                   #                        #           #                #      ",
    "                   #                        #           #                #      ",
    "                   #                        #           #                #      ",
    "                   #########################█###########█#################      ",
    "                                            #           #                       ",
    "                                            #           #                       ",
    "                                            #           #                       ",
    "                                            #           #                       ",
    "                                            #           #                       ",
    "                                            #           #                       ",
    "                                            #############                       ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                   ########     ",
    "                                                                   #      #     ",
    "                                                                   #      #     ",
    "                                                                   #      #     ",
    "                                                                   #      #     ",
    "                                                                   #      #     ",
    "                                                                   ######## ### ",
    "           ############################                                     # # ",
    "           #                          #                                     ### ",
    "           ############################                                         ",
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
]

# Function to count rectangles
def count_rectangles(grid):
    rectangles = 0
    visited = set()
    
    for r in range(len(grid)):
        for c in range(len(grid[r])):
            if grid[r][c] in '#█' and (r, c) not in visited:
                # Start a new rectangle
                rectangles += 1
                # Find the bounds of the rectangle
                end_r, end_c = r, c
                while end_r < len(grid) and grid[end_r][c] in '#█':
                    end_r += 1
                while end_c < len(grid[r]) and grid[r][end_c] in '#█':
                    end_c += 1
                # Mark all parts of the rectangle as visited
                for i in range(r, end_r):
                    for j in range(c, end_c):
                        visited.add((i, j))
    
    return rectangles

# Count the rectangles
total_rectangles = count_rectangles(grid)
print(total_rectangles)
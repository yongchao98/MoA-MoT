def count_rectangles(grid):
    # Determine the maximum length of any line in the grid
    max_length = max(len(line) for line in grid)
    
    # Pad each line to the maximum length
    padded_grid = [line.ljust(max_length) for line in grid]
    
    rows = len(padded_grid)
    cols = len(padded_grid[0]) if rows > 0 else 0
    rectangles = set()

    for r in range(rows):
        for c in range(cols):
            if padded_grid[r][c] in ['#', '█']:
                # Find the bottom-right corner of the rectangle
                br_r, br_c = r, c
                while br_r + 1 < rows and padded_grid[br_r + 1][c] in ['#', '█']:
                    br_r += 1
                while br_c + 1 < cols and padded_grid[r][br_c + 1] in ['#', '█']:
                    br_c += 1
                
                # Check if this is a valid rectangle
                if all(padded_grid[br_r][cc] in ['#', '█'] for cc in range(c, br_c + 1)) and \
                   all(padded_grid[rr][br_c] in ['#', '█'] for rr in range(r, br_r + 1)):
                    rectangles.add((r, c, br_r, br_c))
    
    return len(rectangles)

# Define the grid
grid = [
    "###############                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#             #                                                               ",
    "#          ###█#####################################################         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          #  #                                                    #         ",
    "#          ###█##█████████##############################█████#######         ",
    "#             #  #       █##################            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #   #                ",
    "#             #  #       █                 #            #####                ",
    "#             #  #       █                 #                                 ",
    "#             #  #       █                 #                                 ",
    "#             #  #       █                 #                                 ",
    "###############  #       █                 #                                 ",
    "                    #       █                 #                                 ",
    "                    ########█                 #                                 ",
    "                            #                 #                                 ",
    "                            #                 #                                 ",
    "                           #█#############    #                                 ",
    "                           ##            #    #                                 ",
    "                           ##############█#####                                 ",
    "                           #             #                                      ",
    "                           #             #                                      ",
    "                           #             #       #############################  ",
    "                           #             #       #                           #  ",
    "                           #             #       #                           #  ",
    "                           #             #       #                           #  ",
    "                           #             #       #                           #  ",
    "                           ###############       #                           #  ",
    "                               ###############   #                           #  ",
    "                               #             #   #                           #  ",
    "                               #             #   ███████████████##############  ",
    "                               #             #   #             #                ",
    "                               ###############   ###############                "
]

# Count the rectangles
print(count_rectangles(grid))
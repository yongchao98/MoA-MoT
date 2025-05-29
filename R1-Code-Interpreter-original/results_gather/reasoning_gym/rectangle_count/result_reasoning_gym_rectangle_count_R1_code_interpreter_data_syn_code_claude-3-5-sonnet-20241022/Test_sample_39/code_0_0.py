def is_valid_rectangle(grid, top, left, bottom, right):
    # Check top and bottom edges
    for x in range(left, right + 1):
        if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
            return False
    
    # Check left and right edges
    for y in range(top, bottom + 1):
        if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
            return False
    
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if grid[top][left] not in ['#', '█']:
                continue
                
            for bottom in range(top, height):
                for right in range(left, width):
                    if grid[bottom][right] not in ['#', '█']:
                        continue
                    
                    if is_valid_rectangle(grid, top, left, bottom, right):
                        rectangles.append((top, left, bottom, right))
    
    return len(rectangles)

# Input grid
grid = [
"                                                                                ",
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
"                                                                                "
]

print(find_rectangles(grid))
def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Function to check if a point is part of a rectangle corner
    def is_corner(y, x):
        return (grid[y][x] in ['#', '█']) and \
               (x == 0 or grid[y][x-1] not in ['#', '█']) and \
               (y == 0 or grid[y-1][x] not in ['#', '█'])
    
    # Function to validate rectangle
    def validate_rectangle(top, left, bottom, right):
        # Check all corners exist
        if not all(grid[y][x] in ['#', '█'] for y, x in [
            (top, left), (top, right), (bottom, left), (bottom, right)
        ]):
            return False
        
        # Check all edges are continuous
        for x in range(left, right + 1):
            if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return False
        for y in range(top, bottom + 1):
            if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return False
        return True

    # Store found rectangles to avoid duplicates
    found = set()
    
    # Find all rectangles starting from top-left corners
    for y in range(height):
        for x in range(width):
            if not is_corner(y, x):
                continue
                
            # Find possible right edges
            right = x
            while right < width and grid[y][right] in ['#', '█']:
                right += 1
            
            # Find possible bottom edges
            for r in range(x, right):
                bottom = y
                while bottom < height and grid[bottom][x] in ['#', '█']:
                    bottom += 1
                
                # Check all possible rectangles with these boundaries
                for b in range(y + 1, bottom):
                    for r2 in range(r, right):
                        if validate_rectangle(y, x, b, r2):
                            found.add((y, x, b, r2))

    return len(found)

# Create the grid
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
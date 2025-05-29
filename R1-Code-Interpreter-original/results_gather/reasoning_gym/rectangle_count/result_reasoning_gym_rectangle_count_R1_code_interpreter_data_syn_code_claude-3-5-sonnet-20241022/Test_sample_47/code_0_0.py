# Create the grid as a list of strings
grid = [
"                                                                                ",
"                                                                                ",
"                                                                                ",
"        ###############################################                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        #                                             #                         ",
"        ###############################################                         ",
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
"                              ##################                                ",
"                              #                #                                ",
"                              #                #                                ",
"                              #                #                                ",
"                              ##################                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                        ##############          ",
"                                                        #            #          ",
"                                                        ##############          ",
"                                                                                "
]

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Function to check if a point is part of a rectangle we've already counted
    def is_top_left_corner(i, j):
        if grid[i][j] != '#':
            return False
        # Check if it's not part of an already counted rectangle
        return (i == 0 or grid[i-1][j] != '#') and (j == 0 or grid[i][j-1] != '#')
    
    # For each potential top-left corner
    for i in range(height):
        for j in range(width):
            if is_top_left_corner(i, j):
                # Found a potential rectangle, verify it
                width_count = 0
                height_count = 0
                
                # Count width
                k = j
                while k < width and grid[i][k] == '#':
                    width_count += 1
                    k += 1
                
                # Count height
                k = i
                while k < height and grid[k][j] == '#':
                    height_count += 1
                    k += 1
                
                # Verify it's a complete rectangle
                is_rectangle = True
                for h in range(i, i + height_count):
                    for w in range(j, j + width_count):
                        if h == i or h == i + height_count - 1 or w == j or w == j + width_count - 1:
                            if grid[h][w] != '#':
                                is_rectangle = False
                                break
                        else:
                            if grid[h][w] != ' ':
                                is_rectangle = False
                                break
                    if not is_rectangle:
                        break
                
                if is_rectangle:
                    count += 1

    return count

# Count rectangles and print the result
print(find_rectangles(grid))
def is_complete_rectangle(grid, top, left, bottom, right):
    # Check if all corners are present
    corners = [(top, left), (top, right), (bottom, left), (bottom, right)]
    for y, x in corners:
        if grid[y][x] not in ['#', '█']:
            return False
    
    # Check all edges are continuous
    # Top edge
    for x in range(left + 1, right):
        if grid[top][x] not in ['#', '█']:
            return False
    # Bottom edge
    for x in range(left + 1, right):
        if grid[bottom][x] not in ['#', '█']:
            return False
    # Left edge
    for y in range(top + 1, bottom):
        if grid[y][left] not in ['#', '█']:
            return False
    # Right edge
    for y in range(top + 1, bottom):
        if grid[y][right] not in ['#', '█']:
            return False
    
    # Check if there's empty space inside
    has_space = False
    for y in range(top + 1, bottom):
        for x in range(left + 1, right):
            if grid[y][x] == ' ':
                has_space = True
                break
        if has_space:
            break
    
    return has_space

def find_minimal_rectangles(grid):
    # Convert input to 2D array and remove empty lines
    grid = [list(row) for row in grid.split('\n') if row.strip()]
    
    height = len(grid)
    if height == 0:
        return 0
    width = len(grid[0])
    rectangles = set()
    
    # Find all possible minimal rectangles
    for top in range(height):
        for left in range(width):
            if grid[top][left] not in ['#', '█']:
                continue
            
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if grid[bottom][right] not in ['#', '█']:
                        continue
                    
                    # Check if this forms a valid rectangle
                    if is_complete_rectangle(grid, top, left, bottom, right):
                        # Check if this is not contained within another rectangle
                        is_minimal = True
                        for other_rect in rectangles:
                            t, l, b, r = other_rect
                            if (t <= top and l <= left and b >= bottom and r >= right):
                                is_minimal = False
                                break
                        if is_minimal:
                            rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# Input grid
grid = """                                   ###################################          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
                                   #                                 #          
  #####                            #                                 #          
  #   #                            #                                 #          
  #   #                            #                                 #          
  #   #                            #                                 #          
  #   #                            #                                 #          
  #   #                            #                                 #          
  #   #                            #                                 #          
  #   #                            #                                 #######    
  #   #                            █##############################   ##    #    
  #   #                            █                             #   ##    #    
  #   #                            █                             #   ##    #    
  #   #                            █                             #   ##    #    
  #   #                            █                             #   ##    #    
  #   #                            █                             #   ##    #    
  #   #                            █   ##############            #   ##    #    
  #   #                            █   #            #            #   ##    #    
  #  #█#################           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            #   ##    #    
  #  ##                #           █   #            #            █#####    #    
  #  ##                #           █   #            #            █  ###    #    
  #  ##                #           █   #            #            █  ###    #    
  #  ##                #           █###█############█############█  ###    #    
  #  ##                #           #   #            #            #  ###    #    
  #  ##                #           #   #            #            #  ###    #    
  #  ##                #           #   #            #            #  ###    #    
  #  ##                #           #   #            #            #  ###    #    
  #  ##                #           #   #            #            #  ###    #    
  #  ##                #           #   #            #            ######    #    
  #  ##                #           #   #            #                ##    #    
  #  ##                #           #   ##############                ##    #    
  #  ##                #           #                                 ##    #    
  #  ##            ####█###########█########################         ##    #    
  #  ##            #   #           #                       #         ##    #    
  #  ##            #   #          #█#################      #         ##    #    
  #  ##            #   #          ##                #      #         ##    #    
  #  ##            #   #          #█#################      #         ##    #    
  #  ##        ####█###█###########█##########             #         #######    
  #  #█########█###█####           #         #             #         #          
  #   #        #   #               #         #             #         #          
  #   #        #   #               #         #             #         #          
  #   #        #   ################█#########█##############         #          
  #   #        #                   #         #                       #          
  #   #        #                   #         #                       #          
  #####        #                   #         #                       #          
               #                   #         #                       #          
               #                   ##########█########################          
               #                             #                                  
               #                             #                                  
               #                             #                                  
               #                             #                                  
               #                             #                                  
               #                             #                                  
               #                             #                                  
               #                             #                                  
               ###############################                                  """

print(find_minimal_rectangles(grid))
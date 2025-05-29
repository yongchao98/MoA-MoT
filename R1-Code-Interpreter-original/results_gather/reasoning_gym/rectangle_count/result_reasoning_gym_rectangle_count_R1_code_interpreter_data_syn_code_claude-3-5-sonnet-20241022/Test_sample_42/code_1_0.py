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
    # Convert input to 2D array and remove empty lines
    grid = [list(row) for row in grid.split('\n') if row.strip()]
    
    height = len(grid)
    if height == 0:
        return 0
    width = len(grid[0])
    rectangles = set()
    
    # Find all possible corners (top-left and bottom-right pairs)
    for top in range(height):
        for left in range(width):
            # Skip if not a corner character
            if grid[top][left] not in ['#', '█']:
                continue
                
            # Look for matching bottom-right corners
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    # Skip if not a corner character
                    if grid[bottom][right] not in ['#', '█']:
                        continue
                    
                    # Validate if it forms a complete rectangle
                    if is_valid_rectangle(grid, top, left, bottom, right):
                        # Check if the rectangle has content inside
                        has_inside = False
                        for y in range(top + 1, bottom):
                            for x in range(left + 1, right):
                                if grid[y][x] == ' ':
                                    has_inside = True
                                    break
                            if has_inside:
                                break
                        if has_inside:
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

print(find_rectangles(grid))
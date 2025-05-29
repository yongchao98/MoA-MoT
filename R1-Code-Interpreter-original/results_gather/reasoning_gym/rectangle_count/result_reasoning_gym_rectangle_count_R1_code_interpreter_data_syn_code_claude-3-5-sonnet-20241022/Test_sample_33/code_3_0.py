def is_valid_path(grid, start_x, start_y, end_x, end_y, direction):
    if direction == 'horizontal':
        step = 1 if end_x > start_x else -1
        for x in range(start_x, end_x + step, step):
            if grid[start_y][x] not in '#█':
                return False
    else:  # vertical
        step = 1 if end_y > start_y else -1
        for y in range(start_y, end_y + step, step):
            if grid[y][start_x] not in '#█':
                return False
    return True

def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line]
    height = len(lines)
    width = len(lines[0]) if height > 0 else 0
    
    # Convert to 2D list for easier processing
    grid = [list(line) for line in lines]
    
    rectangle_count = 0
    
    # Find all corner points (positions with '#' or '█')
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                corners.append((x, y))
    
    # For each potential top-left corner
    for top_left in corners:
        x1, y1 = top_left
        
        # Look for potential bottom-right corners
        for bottom_right in corners:
            x2, y2 = bottom_right
            
            # Skip if not a valid rectangle configuration
            if x2 <= x1 or y2 <= y1:
                continue
                
            # Check if we have the other two corners
            if grid[y1][x2] not in '#█' or grid[y2][x1] not in '#█':
                continue
                
            # Verify all four sides
            if (is_valid_path(grid, x1, y1, x2, y1, 'horizontal') and  # top
                is_valid_path(grid, x2, y1, x2, y2, 'vertical') and    # right
                is_valid_path(grid, x1, y2, x2, y2, 'horizontal') and  # bottom
                is_valid_path(grid, x1, y1, x1, y2, 'vertical')):      # left
                
                # Check if interior is empty (except for overlapping parts)
                is_valid = True
                for y in range(y1 + 1, y2):
                    for x in range(x1 + 1, x2):
                        if grid[y][x] == '#':
                            is_valid = False
                            break
                    if not is_valid:
                        break
                
                if is_valid:
                    rectangle_count += 1
    
    return rectangle_count

# The ASCII grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                    ####################################                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #    ####                
         ###########█##################################█####█##█#######         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    ####      #         
         #          #                                  #              #         
         #          #                                  #              #         
         #     #####█#####                             #              #         
         #     #    #    #                             #              #         
         #     #    #    #                             #              #         
         ######█####█####█#############################█###############         
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #              ########  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              ########  
               #    #    #          ###################█#                       
               #    #    #          #                  ##                       
               #    #    #          #                  #█################       
               #    #    #          #                  #█               #       
               #####█#####          #                  #█               #       
                    ################█████████████████████               #       
                                                        #               #       
                                                        #               #       
                                                        #               #       
                                                        #               #       
                                                        #               #       
                                                        #################       
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))
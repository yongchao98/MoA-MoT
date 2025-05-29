def count_rectangles(grid):
    lines = grid.strip('\n').split('\n')
    height = len(lines)
    width = max(len(line) for line in lines)
    
    def is_valid_edge_char(i, j):
        if i < 0 or i >= height or j < 0 or j >= width or j >= len(lines[i]):
            return False
        return lines[i][j] in '#█'
    
    def find_rectangle_from_corner(start_i, start_j):
        # Must be a corner character
        if not is_valid_edge_char(start_i, start_j):
            return None
            
        # Must not have edge character above or to the left (to ensure top-left corner)
        if is_valid_edge_char(start_i-1, start_j) or is_valid_edge_char(start_i, start_j-1):
            return None
            
        # Find right edge
        right = start_j + 1
        while right < width and right < len(lines[start_i]) and is_valid_edge_char(start_i, right):
            right += 1
        right -= 1
        
        if right <= start_j:
            return None
            
        # Find bottom edge
        bottom = start_i + 1
        while bottom < height and is_valid_edge_char(bottom, start_j):
            bottom += 1
        bottom -= 1
        
        if bottom <= start_i:
            return None
            
        # Verify bottom-right corner exists
        if not is_valid_edge_char(bottom, right):
            return None
            
        # Verify all edges are complete
        for i in range(start_i, bottom + 1):
            if not is_valid_edge_char(i, start_j) or not is_valid_edge_char(i, right):
                return None
                
        for j in range(start_j, right + 1):
            if not is_valid_edge_char(start_i, j) or not is_valid_edge_char(bottom, j):
                return None
                
        return (start_i, start_j, bottom, right)
    
    rectangles = set()
    # Only look for top-left corners
    for i in range(height):
        for j in range(len(lines[i])):
            rect = find_rectangle_from_corner(i, j)
            if rect:
                rectangles.add(rect)
    
    # Debug: Print rectangle coordinates
    for rect in rectangles:
        top, left, bottom, right = rect
        print(f"Rectangle: ({top}, {left}) to ({bottom}, {right})")
    
    print(len(rectangles))

# Your ASCII grid
grid = """    #################################################################           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                              #################█########   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                   ###########████#############█########   
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             ########    
    ####################################█#############█###############     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        ###############              #     #    
                                                                     #######    
                                                                      ########  
                                                                      #      #  
                                                                      ########  """

count_rectangles(grid)
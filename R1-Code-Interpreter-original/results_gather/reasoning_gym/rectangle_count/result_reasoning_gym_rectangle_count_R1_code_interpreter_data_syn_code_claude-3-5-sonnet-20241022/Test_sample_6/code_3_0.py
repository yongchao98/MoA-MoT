def find_rectangles(grid):
    # Convert grid to list of lines
    lines = grid.split('\n')
    rows = len(lines)
    cols = len(lines[0]) if rows > 0 else 0
    
    def is_border(char):
        return char in ['#', '█']
    
    def validate_rectangle(top, left, height, width):
        # Check top and bottom borders
        top_border = lines[top][left:left+width]
        bottom_border = lines[top+height-1][left:left+width]
        if not all(is_border(c) for c in top_border) or not all(is_border(c) for c in bottom_border):
            return False
        
        # Check left and right borders
        for y in range(top, top+height):
            if not is_border(lines[y][left]) or not is_border(lines[y][left+width-1]):
                return False
            
        # Check if there's content inside
        for y in range(top+1, top+height-1):
            for x in range(left+1, left+width-1):
                if lines[y][x] not in [' ', '#', '█']:
                    return False
        
        return True
    
    rectangles = []
    # Find rectangles by scanning for top-left corners
    for y in range(rows-1):
        for x in range(cols-1):
            if is_border(lines[y][x]):
                # Look for matching bottom-right corners
                for h in range(2, rows-y+1):
                    for w in range(2, cols-x+1):
                        if y+h <= rows and x+w <= cols:
                            if validate_rectangle(y, x, h, w):
                                # Check if this is a minimal rectangle
                                is_minimal = True
                                for inner_y in range(y+1, y+h-1):
                                    for inner_x in range(x+1, x+w-1):
                                        if is_border(lines[inner_y][inner_x]):
                                            is_minimal = False
                                            break
                                    if not is_minimal:
                                        break
                                if is_minimal:
                                    rectangles.append((y, x, h, w))
    
    # Count distinct rectangles
    distinct_rectangles = set()
    for rect in rectangles:
        y, x, h, w = rect
        # Create a signature for the rectangle
        signature = tuple((y+i, x+j) for i in range(h) for j in range(w) 
                         if is_border(lines[y+i][x+j]))
        distinct_rectangles.add(signature)
    
    return len(distinct_rectangles)

# Input grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                      ##########################
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      ##########################
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                       ###################      
                                                       #                 #      
                                                       ###################      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                         ####################                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #    ##############█#################  
                                         #    #             #                #  
                                         #    #             #                #  
                                         #    ##############█#################  
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
          ###############################█##################█###################
          #                              ####################                  #
          #                                                                    #
          #                                                                    #
          #                                                                    #
          #                                                                    #
          ######################################██████████████##################
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #               ###
                                                #            #               # #
                                                #            #               # #
                                                #            #               ###
                                                ##############                  """

print(find_rectangles(grid))
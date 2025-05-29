def find_rectangles(grid):
    # Convert grid to list of strings and clean empty lines
    lines = [line for line in grid.split('\n') if line]
    
    def is_valid_char(char):
        return char in '#█'
    
    def trace_rectangle(start_y, start_x):
        # Find the width of top edge
        width = 0
        x = start_x
        while x < len(lines[start_y]) and is_valid_char(lines[start_y][x]):
            width += 1
            x += 1
            
        if width < 2:  # Too narrow to be a rectangle
            return False
            
        # Find the height of left edge
        height = 0
        y = start_y
        while y < len(lines) and is_valid_char(lines[y][start_x]):
            height += 1
            y += 1
            
        if height < 2:  # Too short to be a rectangle
            return False
            
        # Verify bottom edge
        x = start_x
        for i in range(width):
            if not (start_x + i < len(lines[start_y + height - 1]) and 
                   is_valid_char(lines[start_y + height - 1][start_x + i])):
                return False
                
        # Verify right edge
        y = start_y
        for i in range(height):
            if not (start_x + width - 1 < len(lines[start_y + i]) and 
                   is_valid_char(lines[start_y + i][start_x + width - 1])):
                return False
                
        # Verify there's empty space inside
        has_space = False
        for y in range(start_y + 1, start_y + height - 1):
            for x in range(start_x + 1, start_x + width - 1):
                if lines[y][x] == ' ':
                    has_space = True
                    break
            if has_space:
                break
                
        return has_space
    
    count = 0
    # Scan for potential top-left corners
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            if is_valid_char(lines[y][x]):
                # Check if this could be a top-left corner
                if (x == 0 or not is_valid_char(lines[y][x-1])) and \
                   (y == 0 or not is_valid_char(lines[y-1][x])):
                    if trace_rectangle(y, x):
                        count += 1
    
    return count

# Your grid
grid = """                                                                                
                                                                                
                                                                                
                                 ###########################                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
   ##############################█#########                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #########█#################                    
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                            ##########█###################################  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #############################█##########                                  #  
                                ##############################################  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                            ############################                        
                            #                          #                        
                            ############################                        
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))
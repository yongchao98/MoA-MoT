def count_rectangles(grid):
    # Convert grid to list of strings
    lines = grid.split('\n')
    # Remove empty lines and ensure all lines have same width
    lines = [line for line in lines if line.strip()]
    
    def is_valid_char(char):
        return char in '#█'
    
    def verify_rectangle(top, left, bottom, right):
        # Check if all corners exist
        if not is_valid_char(lines[top][left]) or \
           not is_valid_char(lines[top][right]) or \
           not is_valid_char(lines[bottom][left]) or \
           not is_valid_char(lines[bottom][right]):
            return False
            
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_valid_char(lines[top][x]) or \
               not is_valid_char(lines[bottom][x]):
                return False
                
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_valid_char(lines[y][left]) or \
               not is_valid_char(lines[y][right]):
                return False
                
        # Check if interior is consistent
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if is_valid_char(lines[y][x]) and lines[y][x] != ' ':
                    return False
        
        return True
    
    # Find distinct rectangles
    rectangles = set()
    height = len(lines)
    width = len(lines[0])
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if is_valid_char(lines[top][left]):
                for bottom in range(top + 1, height):
                    for right in range(left + 1, width):
                        if verify_rectangle(top, left, bottom, right):
                            rectangles.add((top, left, bottom, right))
    
    # Count base rectangles
    base_count = len(rectangles)
    
    # Count overlaps (each █ represents one additional rectangle)
    overlap_count = sum(1 for line in lines for char in line if char == '█')
    
    total = base_count + overlap_count
    print(total)

grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
          ####################################################                  
          #                                                  #                  
          #                                                  #                  
          #              ########                            #                  
          #              #      #                            #                  
          #            ##█######█#####                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # ########    #                       #                  
          #            ###############                       #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          ####################################################                  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

count_rectangles(grid)
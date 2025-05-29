def count_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    
    def is_border(char):
        return char in '#█'
    
    def find_complete_rectangle(y, x):
        # Find width of rectangle
        width = 0
        for w in range(x, len(lines[y])):
            if not is_border(lines[y][w]):
                break
            width += 1
        if width < 2:  # Need at least 2 characters for width
            return None
            
        # Find height of rectangle
        height = 0
        for h in range(y, len(lines)):
            if not is_border(lines[h][x]):
                break
            # Check right border
            if not is_border(lines[h][x + width - 1]):
                return None
            height += 1
        if height < 2:  # Need at least 2 characters for height
            return None
            
        # Verify top and bottom borders are complete
        for w in range(x, x + width):
            if not is_border(lines[y][w]) or not is_border(lines[y + height - 1][w]):
                return None
                
        # Verify left and right borders are complete
        for h in range(y, y + height):
            if not is_border(lines[h][x]) or not is_border(lines[h][x + width - 1]):
                return None
                
        return (y, x, height, width)
    
    rectangles = set()
    overlap_count = 0
    
    # Find all rectangles
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            if is_border(lines[y][x]):
                rect = find_complete_rectangle(y, x)
                if rect:
                    rectangles.add(rect)
    
    # Count overlaps (█ characters)
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            if lines[y][x] == '█':
                overlap_count += 1
    
    # Each █ represents two overlapping rectangles, but we only need to add one
    # since the base rectangles are already counted
    total = len(rectangles) + overlap_count
    print(total)

# Create the grid string (your ASCII grid here)
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
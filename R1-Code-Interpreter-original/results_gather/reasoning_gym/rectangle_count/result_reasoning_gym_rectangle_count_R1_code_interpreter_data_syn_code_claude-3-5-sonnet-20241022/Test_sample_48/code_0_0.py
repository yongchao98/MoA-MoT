def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = max(len(line) for line in lines)
    
    # Standardize width of all lines
    lines = [line.ljust(width) for line in lines]
    
    def is_corner(i, j):
        if lines[i][j] not in '#█':
            return False
        # Check if it's a corner by looking at adjacent cells
        if i > 0 and i < height - 1 and j > 0 and j < width - 1:
            return ((lines[i][j-1] in '#█' or lines[i][j+1] in '#█') and 
                    (lines[i-1][j] in '#█' or lines[i+1][j] in '#█'))
        return False
    
    def trace_rectangle(top, left):
        # Find bottom right corner of rectangle
        right = left
        bottom = top
        
        # Find right edge
        while right + 1 < width and lines[top][right + 1] in '#█':
            right += 1
            
        # Find bottom edge
        while bottom + 1 < height and lines[bottom + 1][left] in '#█':
            bottom += 1
            
        # Verify rectangle
        for i in range(top, bottom + 1):
            if lines[i][left] not in '#█' or lines[i][right] not in '#█':
                return None
        for j in range(left, right + 1):
            if lines[top][j] not in '#█' or lines[bottom][j] not in '#█':
                return None
                
        return (top, left, bottom, right)
    
    rectangles = set()
    # Find all rectangles by looking for top-left corners
    for i in range(height):
        for j in range(width):
            if is_corner(i, j):
                rect = trace_rectangle(i, j)
                if rect:
                    rectangles.add(rect)
    
    # Count overlapping rectangles (marked with █)
    overlap_count = sum(1 for i in range(height) 
                       for j in range(width) 
                       if lines[i][j] == '█')
    
    return len(rectangles) + overlap_count

# Test with the provided grid
grid = """                                                                                
                                                                   ############ 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   #          # 
                                                                   ############ 
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                      ######### 
                                                                      #       # 
      ####################################################            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      #                                                  #            #       # 
      ####################################################            #       # 
                                                                      #       # 
                                                                      #       # 
                                                                      #       # 
                                                                      #       # 
                                              ########################█###    # 
                                              #                       #  #    # 
                                              ########################█###    # 
                                                                      #       # 
                                                                      #       # 
                                                                      #       # 
                                                  #####################       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #                  ##       # 
                                                  #####################       # 
                                                                      #       # 
                                                                      #       # 
                                                                      ######### 
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))
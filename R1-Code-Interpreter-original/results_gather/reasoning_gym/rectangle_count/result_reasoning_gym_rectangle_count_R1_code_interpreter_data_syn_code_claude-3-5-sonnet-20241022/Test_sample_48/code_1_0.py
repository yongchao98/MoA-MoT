def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = max(len(line) for line in lines)
    
    # Standardize width of all lines
    lines = [line.ljust(width) for line in lines]
    
    def is_top_left_corner(i, j):
        if lines[i][j] not in '#█':
            return False
        # Check if it's a top-left corner
        left_clear = j == 0 or lines[i][j-1] not in '#█'
        top_clear = i == 0 or lines[i-1][j] not in '#█'
        return left_clear and top_clear and (j+1 < width and lines[i][j+1] in '#█') and (i+1 < height and lines[i+1][j] in '#█')
    
    def find_rectangle(top, left):
        right = left
        bottom = top
        
        # Find right edge
        while right + 1 < width and lines[top][right + 1] in '#█':
            right += 1
            
        # Find bottom edge
        while bottom + 1 < height and lines[bottom + 1][left] in '#█':
            bottom += 1
            
        # Verify it's a complete rectangle
        for i in range(top, bottom + 1):
            if lines[i][left] not in '#█' or lines[i][right] not in '#█':
                return None
        for j in range(left, right + 1):
            if lines[top][j] not in '#█' or lines[bottom][j] not in '#█':
                return None
                
        return (top, left, bottom, right)
    
    rectangles = set()
    # Find rectangles by looking for top-left corners
    for i in range(height):
        for j in range(width):
            if is_top_left_corner(i, j):
                rect = find_rectangle(i, j)
                if rect:
                    rectangles.add(rect)
    
    # Count overlapping areas (marked with █)
    overlap_count = 0
    for i in range(height):
        for j in range(width):
            if lines[i][j] == '█':
                # Check if this █ is part of an overlap between two rectangles
                top = bottom = i
                left = right = j
                # Expand to find the full overlapping area
                while top > 0 and lines[top-1][j] == '█':
                    top -= 1
                while bottom < height-1 and lines[bottom+1][j] == '█':
                    bottom += 1
                while left > 0 and lines[i][left-1] == '█':
                    left -= 1
                while right < width-1 and lines[i][right+1] == '█':
                    right += 1
                if (top, left) == (i, j):  # Only count each overlapping area once
                    overlap_count += 1
    
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
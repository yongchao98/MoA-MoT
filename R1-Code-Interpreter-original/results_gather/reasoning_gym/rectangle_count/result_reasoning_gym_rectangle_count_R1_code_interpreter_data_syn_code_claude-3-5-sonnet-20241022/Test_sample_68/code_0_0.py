def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Function to check if a point is part of rectangle border
    def is_border(ch):
        return ch in '#█'
    
    # Find all rectangles
    for y in range(height):
        for x in range(width):
            # If we find a top-left corner
            if is_border(grid[y][x]):
                # Skip if this point is already part of a found rectangle
                if any(r[0] <= x <= r[2] and r[1] <= y <= r[3] for r in rectangles):
                    continue
                
                # Find right edge
                right = x
                while right < width and is_border(grid[y][right]):
                    right += 1
                right -= 1
                
                # Find bottom edge
                bottom = y
                while bottom < height and is_border(grid[bottom][x]):
                    bottom += 1
                bottom -= 1
                
                # Verify rectangle
                valid = True
                for cy in range(y, bottom + 1):
                    if not (is_border(grid[cy][x]) and is_border(grid[cy][right])):
                        valid = False
                        break
                for cx in range(x, right + 1):
                    if not (is_border(grid[y][cx]) and is_border(grid[bottom][cx])):
                        valid = False
                        break
                
                if valid:
                    rectangles.append((x, y, right, bottom))
    
    # Count overlapping points to determine if we need to add extra rectangles
    overlap_count = 0
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                overlap_count += 1
    
    # Each '█' character represents an overlap of two rectangles
    total_rectangles = len(rectangles) + overlap_count
    
    print(total_rectangles)

# Create the grid string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                ############### 
                                                                #             # 
                                                                #             # 
                                                                #             # 
                                                                #             # 
                                                                #             # 
                                                                #             # 
                                                                #             # 
                                                          ######█##########   # 
                                                          #     #         #   # 
                                                          #     #         #   # 
                                                          #     #         #   # 
                                                          #     #         #   # 
                                                          #     #         #   # 
                                                          #     #         #   # 
                                                          #     #         #   # 
                                                          #     #         #   # 
                                                          #     ##########█#### 
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                          #               #     
                                                  ########█###############█#    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       #               ##    
                                                  #       ##################    
                                                  #                        #    
                                                  #                        #    
                                                  #                        #    
                                                  ##########################    
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)
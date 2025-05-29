def find_rectangles(grid):
    # Convert grid to list of strings
    lines = grid.split('\n')
    # Remove empty lines
    grid = [line for line in lines if line.strip()]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    # Find all corners (top-left corners)
    corners = []
    for i in range(height-1):
        for j in range(width-1):
            # Check if this is a top-left corner of a rectangle
            if grid[i][j] in '#█':
                # Look right for the top-right corner
                right = j + 1
                while right < width and grid[i][right] in '#█':
                    right += 1
                right -= 1
                
                # Look down for the bottom-left corner
                bottom = i + 1
                while bottom < height and grid[bottom][j] in '#█':
                    bottom += 1
                bottom -= 1
                
                # Verify it's a complete rectangle
                is_rectangle = True
                # Check bottom edge
                for x in range(j, right + 1):
                    if grid[bottom][x] not in '#█':
                        is_rectangle = False
                        break
                # Check right edge
                for y in range(i, bottom + 1):
                    if grid[y][right] not in '#█':
                        is_rectangle = False
                        break
                
                if is_rectangle:
                    corners.append((i, j, bottom, right))
    
    # Count overlapping points
    overlap_count = 0
    for i in range(height):
        for j in range(width):
            if grid[i][j] == '█':
                overlap_count += 1
    
    # Each overlap point represents where two rectangles were counted as one
    total_rectangles = len(corners) - overlap_count + overlap_count
    print(total_rectangles)

# Create test grid
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
                                                 #          #             ###   
                                                 #          #             # #   
                                                 #          #             # #   
                                             ####█########  #             # #   
                                             #   #       #  #             # #   
                                             #   #       #  #             # #   
                                             #   #       #  #             # #   
                                             ####█████████###             # #   
                                                                          # #   
                                                                          # #   
                                                                          ###   
                                                                                
                                                                                """

find_rectangles(grid)
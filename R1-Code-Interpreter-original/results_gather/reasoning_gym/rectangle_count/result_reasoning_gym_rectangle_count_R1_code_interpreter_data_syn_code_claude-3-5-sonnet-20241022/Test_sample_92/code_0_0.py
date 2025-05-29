def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Helper function to check if a point is a corner
    def is_corner(y, x):
        if y < 0 or y >= height or x < 0 or x >= width:
            return False
        return grid[y][x] in ['#', '█']
    
    # Helper function to validate rectangle
    def validate_rectangle(top, left, bottom, right):
        # Check corners
        if not all(is_corner(y, x) for y, x in [
            (top, left), (top, right),
            (bottom, left), (bottom, right)
        ]):
            return False
            
        # Check horizontal edges
        for x in range(left + 1, right):
            if not (grid[top][x] in ['#', '█'] and grid[bottom][x] in ['#', '█']):
                return False
                
        # Check vertical edges
        for y in range(top + 1, bottom):
            if not (grid[y][left] in ['#', '█'] and grid[y][right] in ['#', '█']):
                return False
                
        return True

    # Convert input string to list of lines
    grid = [list(line) for line in grid.strip().split('\n')]
    
    # Find all potential corners
    corners = [(y, x) for y in range(height) for x in range(width) 
              if grid[y][x] in ['#', '█']]
    
    # Try all possible rectangle combinations
    for i, (top, left) in enumerate(corners):
        for bottom, right in corners[i:]:
            if bottom > top and right > left:  # Ensure correct orientation
                if validate_rectangle(top, left, bottom, right):
                    rectangles.append((top, left, bottom, right))
    
    return len(rectangles)

# Test grid (your provided grid)
test_grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                      ###       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      ###       
                                                                                
                                                      ###                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      ###                       
                                                                                
                                                     ###########                
                                                     #         #                
                                                     #         #                
                                                     #         #                
                                        ########################                
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #############                           
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                     ###                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     ###                        
                                                                                """

print(find_rectangles(test_grid))
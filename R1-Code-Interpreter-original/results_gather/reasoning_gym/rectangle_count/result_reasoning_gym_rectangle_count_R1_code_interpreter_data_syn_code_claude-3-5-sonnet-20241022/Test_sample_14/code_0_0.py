def find_rectangles(grid):
    # Convert string grid to list of lists for easier processing
    grid = [list(row) for row in grid.splitlines()]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    rectangles = 0
    
    # Function to check if point is part of rectangle border
    def is_border(i, j):
        return 0 <= i < height and 0 <= j < width and (grid[i][j] == '#' or grid[i][j] == '█')
    
    # For each point in grid
    for i in range(height):
        for j in range(width):
            # If we find a top-left corner
            if is_border(i, j) and not is_border(i-1, j) and not is_border(i, j-1):
                # Look for width of rectangle
                width_rect = 0
                while j + width_rect < width and is_border(i, j + width_rect):
                    width_rect += 1
                
                if width_rect < 2:  # Too narrow to be rectangle
                    continue
                    
                # Look for height of rectangle
                height_rect = 0
                while i + height_rect < height and is_border(i + height_rect, j):
                    height_rect += 1
                    
                if height_rect < 2:  # Too short to be rectangle
                    continue
                
                # Verify if it's a complete rectangle
                is_complete = True
                # Check bottom edge
                for w in range(width_rect):
                    if not is_border(i + height_rect - 1, j + w):
                        is_complete = False
                        break
                # Check right edge
                for h in range(height_rect):
                    if not is_border(i + h, j + width_rect - 1):
                        is_complete = False
                        break
                
                if is_complete:
                    rectangles += 1
    
    print(rectangles)

# Create the grid string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                             ###########################################        
                             #                                         #        
                             #                                         #        
                             #                                         #        
                             ###########################################        
                                                                                
                                                                                
                                                                                
                                                                                
                                                   ####                         
                                                   #  #                         
                                                   #  #                         
                                                   #  #                         
                                                   ####                         """

find_rectangles(grid)
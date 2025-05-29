def count_rectangles(grid_str):
    # Split into lines and remove any empty lines at the start/end
    grid = [line for line in grid_str.splitlines() if line]
    height = len(grid)
    width = len(grid[0]) if grid else 0
    
    # Count distinct rectangles and overlaps
    rectangles = 0
    overlaps = 0
    
    # Helper function to validate if a point is part of rectangle border
    def is_border(y, x):
        return 0 <= y < height and 0 <= x < width and (grid[y][x] in '#█')
    
    # Helper function to verify if a point is top-left corner of a rectangle
    def is_top_left_corner(y, x):
        if not is_border(y, x):
            return False
        # Check if there's no border above and to the left
        return (not is_border(y-1, x) and not is_border(y, x-1))
    
    # For each potential top-left corner
    for y in range(height):
        for x in range(width):
            if is_top_left_corner(y, x):
                rectangles += 1
                
            # Count overlaps separately
            if grid[y][x] == '█':
                overlaps += 1
    
    # Each overlap represents an additional rectangle
    total = rectangles + overlaps
    
    # Debug print
    print(f"Found {rectangles} base rectangles and {overlaps} overlaps")
    print(f"Total: {total}")

# Create the exact grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                    ####███#    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    ####█#█#    
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        ###     
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
      ######################                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                      #######       
      #                    #                                      #     #       
      ######################                                      #######       
                                                                                
                                                                                
                                                                                """

count_rectangles(grid)
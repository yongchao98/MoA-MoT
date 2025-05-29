def count_rectangles(grid_str):
    # Split into lines and remove any empty lines at the start/end
    grid = [line for line in grid_str.splitlines() if line]
    
    # Initialize counters
    basic_rectangles = 0
    overlaps = 0
    
    # First, let's count all the basic rectangles by finding their top-left corners
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            # Check if current position is a top-left corner of a rectangle
            if grid[i][j] in '#█':
                # Check if it's a top-left corner by looking at positions above and left
                is_corner = True
                if i > 0 and grid[i-1][j] in '#█':  # Check above
                    is_corner = False
                if j > 0 and grid[i][j-1] in '#█':  # Check left
                    is_corner = False
                
                if is_corner:
                    basic_rectangles += 1
            
            # Count overlaps (each █ represents two rectangles)
            if grid[i][j] == '█':
                overlaps += 1
    
    # Each overlap point (█) means we need to add one more rectangle
    # (as it's already counted once in basic_rectangles)
    total_rectangles = basic_rectangles + overlaps
    
    print(total_rectangles)

# Create the exact grid with the correct characters
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
def count_rectangles(grid):
    # Convert the grid into a list of strings
    grid_lines = grid.strip().split('\n')
    
    # Store the rectangles as tuples of (top-left, bottom-right)
    rectangles = set()
    
    # Iterate over each line to find horizontal lines
    for i, line in enumerate(grid_lines):
        # Find all horizontal lines in the current line
        start = None
        for j, char in enumerate(line):
            if char in '#█':
                if start is None:
                    start = j
            else:
                if start is not None:
                    # We found the end of a horizontal line
                    end = j - 1
                    # Now, find the vertical lines to complete the rectangle
                    for k in range(i + 1, len(grid_lines)):
                        if grid_lines[k][start] in '#█' and grid_lines[k][end] in '#█':
                            # Check if this forms a rectangle
                            if all(grid_lines[k][x] in '#█' for x in range(start, end + 1)):
                                # Add the rectangle to the set
                                rectangles.add((i, start, k, end))
                        else:
                            break
                    start = None
    
    return len(rectangles)

# Define the grid
grid = """
                ###################################                             
                #                          #####  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #####  #                             
                #                       ##########█#################            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       ##########█#################            
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 # #######################     
                #                                 # #                     #     
                #                                 # #                     #     
                #                                 # #                     #     
                #                                 # #                     #     
                #               ##################█#█###                  #     
                #               #                 # #  #                  #     
                #               #                 # ###█###################     
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #     ################## 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               ##################█#####     #                # 
                #                                 #          #                # 
                #                                 #          #                # 
                #                                 #          #                # 
                #                                 #          ################## 
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                ###################################                             
"""

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
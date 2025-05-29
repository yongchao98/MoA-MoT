def find_rectangles(grid):
    lines = grid.splitlines()
    count = 0
    
    # Find top-left corners of rectangles
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            # Check if this is a top-left corner
            if lines[y][x] in '#█':
                # Check if it's actually a corner (empty space or grid edge to the left and above)
                is_corner = (x == 0 or lines[y][x-1] not in '#█') and (y == 0 or lines[y-1][x] not in '#█')
                
                if is_corner:
                    # If we found a corner, check if it's part of an overlapping section
                    has_overlap = False
                    # Look for █ in the current position
                    if lines[y][x] == '█':
                        has_overlap = True
                    count += 2 if has_overlap else 1
    
    return count

# Test with the provided ASCII art
ascii_art = """                                                                                
                                                                                
                                      ########################                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      ########################                  
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
                                                             # #                
                                                             # #                
                                                             # #                
                                                     ########█#█#######         
                                                     #       # #      #         
        ####################                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         ########█#█#######         
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        ####################                                 # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             ###                
                                                                                
                           ###################################################  
                           #                                                 #  
                           #                                                 #  
                           #                                                 #  
                           ###################################################  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(ascii_art))
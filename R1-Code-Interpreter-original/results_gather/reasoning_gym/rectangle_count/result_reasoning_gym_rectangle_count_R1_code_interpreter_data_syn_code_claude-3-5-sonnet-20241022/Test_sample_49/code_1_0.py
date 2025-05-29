def find_rectangles(grid):
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    def is_border(char):
        return char in '#█'
    
    def find_complete_rectangle(start_y, start_x):
        # Skip if not a potential start of rectangle
        if not is_border(lines[start_y][start_x]):
            return False
            
        # Skip if not likely to be top-left corner
        if start_y > 0 and is_border(lines[start_y-1][start_x]):
            return False
        if start_x > 0 and is_border(lines[start_y][start_x-1]):
            return False
            
        # Find width of top edge
        right = start_x
        while right < width and is_border(lines[start_y][right]):
            right += 1
        right -= 1
        
        if right == start_x:  # Not wide enough to be a rectangle
            return False
            
        # Find height of left edge
        bottom = start_y
        while bottom < height and is_border(lines[bottom][start_x]):
            bottom += 1
        bottom -= 1
        
        if bottom == start_y:  # Not tall enough to be a rectangle
            return False
            
        # Verify bottom and right edges exist
        for x in range(start_x, right + 1):
            if not is_border(lines[bottom][x]):
                return False
                
        for y in range(start_y, bottom + 1):
            if not is_border(lines[y][right]):
                return False
                
        # Check if this rectangle contains any '█' characters
        has_overlap = False
        for y in range(start_y, bottom + 1):
            for x in range(start_x, right + 1):
                if lines[y][x] == '█':
                    has_overlap = True
                    break
            if has_overlap:
                break
                
        # If there's an overlap character, this indicates two rectangles
        return 2 if has_overlap else 1
    
    total_rectangles = 0
    # Scan for all potential top-left corners
    for y in range(height):
        for x in range(width):
            rect_count = find_complete_rectangle(y, x)
            if rect_count:
                total_rectangles += rect_count
                
    # Adjust for double counting of overlapping rectangles
    total_rectangles = total_rectangles // 2
    
    return total_rectangles

# Test grid (same as before)
grid = """
                                                                                
                                                                                
                                                 ##################             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 ##################             
                                                                                
                                                                                
                                                                                
                                                                                
                                                 ###############################
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                             #
                                                 #                     ####    #
                                                 #                     #  #    #
                                                 #                     #  #    #
                                        #########█###########          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          ████#####
                                        #        #          #          #      ##
                                        #        #       ###█######    #      ##
                                        #        #       #  #     #    #      ##
                                        #        #       ###█######    #      ##
                                        #        #          #          #      ##
                                        #########█###########          #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #########
                                                 #                             #
                                                 #                             #
                                                 #        #############        #
                                                 #########█###########█#########
                                                          #           #         
                                                          #           #         
                                                          #           #         
                                                          #           #         
                   ####################███████████        #############         
                   #                   #         █                              
                   ####################███████████                              
"""

print(find_rectangles(grid))
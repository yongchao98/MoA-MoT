def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = len(rows[0])
    
    def is_top_left_corner(r, c):
        return rows[r][c] in '#█' and (r == 0 or rows[r-1][c] == ' ') and (c == 0 or rows[r][c-1] == ' ')
    
    def is_bottom_right_corner(r, c):
        return rows[r][c] in '#█' and (r == height-1 or rows[r+1][c] == ' ') and (c == width-1 or rows[r][c+1] == ' ')
    
    top_left_corners = []
    bottom_right_corners = []
    
    for r in range(height):
        for c in range(width):
            if is_top_left_corner(r, c):
                top_left_corners.append((r, c))
            if is_bottom_right_corner(r, c):
                bottom_right_corners.append((r, c))
    
    rectangle_count = 0
    for (r1, c1) in top_left_corners:
        for (r2, c2) in bottom_right_corners:
            if r1 < r2 and c1 < c2:
                # Check if this forms a valid rectangle
                if all(rows[r1][c] in '#█' for c in range(c1, c2+1)) and \
                   all(rows[r2][c] in '#█' for c in range(c1, c2+1)) and \
                   all(rows[r][c1] in '#█' for r in range(r1, r2+1)) and \
                   all(rows[r][c2] in '#█' for r in range(r1, r2+1)):
                    rectangle_count += 1
    
    return rectangle_count

grid = """
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
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #    #####################     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                #############█####█#######            #     
                   #                #            #    #      #            #     
                   #                #            #    #######█#############     
                   #                #            #           #                  
                   #                #            #           #                  
                   #                #            #           #                  
                   #                #            #           #              ### 
                   #                #            #           #              # # 
                   #                #############█############              # # 
                   #                             #                          # # 
                   #                             #                          # # 
                   #                             #                          # # 
                   #                             #  ####################    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
     ##############█##                           #  #                  #    # # 
     #             # #                           #  #                  #    # # 
     #             # #                           #  #                  #    # # 
     ##############█##                           #  ####################    # # 
                   #               ##############█##########################██# 
                   ################█##############                           #  
                                   #                                         #  
                                   ###########################################  
"""

print(count_rectangles(grid))
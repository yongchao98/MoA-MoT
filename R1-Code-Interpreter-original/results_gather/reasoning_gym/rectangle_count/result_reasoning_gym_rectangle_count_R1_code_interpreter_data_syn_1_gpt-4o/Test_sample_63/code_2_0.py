def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    num_cols = len(rows[0]) if num_rows > 0 else 0

    def is_rectangle_char(c):
        return c == '#' or c == '█'

    rectangles = set()

    # Scan for horizontal lines
    for row in range(num_rows):
        in_rectangle = False
        start_col = 0
        for col in range(num_cols):
            if is_rectangle_char(rows[row][col]):
                if not in_rectangle:
                    in_rectangle = True
                    start_col = col
            else:
                if in_rectangle:
                    in_rectangle = False
                    rectangles.add((row, start_col, row, col - 1))
        if in_rectangle:
            rectangles.add((row, start_col, row, num_cols - 1))

    # Scan for vertical lines
    for col in range(num_cols):
        in_rectangle = False
        start_row = 0
        for row in range(num_rows):
            if is_rectangle_char(rows[row][col]):
                if not in_rectangle:
                    in_rectangle = True
                    start_row = row
            else:
                if in_rectangle:
                    in_rectangle = False
                    rectangles.add((start_row, col, row - 1, col))
        if in_rectangle:
            rectangles.add((start_row, col, num_rows - 1, col))

    # Count each '█' as two rectangles
    overlap_count = sum(row.count('█') for row in rows)
    total_rectangles = len(rectangles) + overlap_count

    return total_rectangles

grid = """
                                                            ############        
                                                            #          #        
                                                            #          #        
                                                            #          #        
                                                            #          #        
                                                            ############        
                                                                                
                                                                                
                                                                                
 ###############################################################################
 #                                                                             #
 #                  ##########                                                 #
 #                  #        #                                                 #
 #                  #        #                                                 #
 #                  #        #                                                 #
 #                  #        #                                                 #
 #                  #        #                                 ###             #
 #                  #        #                                 # #             #
 #                  #        #                                 # #             #
 #                  #        #                                 # #             #
 #                  #        #                                 # #             #
 #                  #        #                 #######         # #             #
 #                  #        #                 #     #         # #             #
 #                  #        #                 #     #         # # #############
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #                  #        #                 #     #         # # #          ##
 #           #######█########█#################███████######## # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 #           #      #        #                               # # # #          ##
 ############█######█########█###############################█#█#█#█##########█#
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      #        #                               # # # #          # 
             #      ##########                               # # # #          # 
             ################################################# # # #          # 
                                                               # # #          # 
                                                               # # #          # 
                                                               # # #          # 
                                                               # # #          # 
                                                               # # ############ 
                                                               # #              
                                                               # #              
                                                               # █############# 
                                                               ##█            # 
                                                                 #            # 
                        #####                                    #            # 
                        #   #                                    ############## 
                        #   #                                                   
                        #   #                                                   
                        #   #                                                   
                        #   #                                                   
                        #   #                                                   
                        #####                                                   
"""

print(count_rectangles(grid))
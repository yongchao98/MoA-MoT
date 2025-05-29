def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    num_cols = len(rows[0]) if num_rows > 0 else 0

    def is_rectangle_char(c):
        return c == '#' or c == '█'

    def find_rectangles_in_row(row):
        rectangles = []
        in_rectangle = False
        start_col = 0
        for col in range(num_cols):
            if is_rectangle_char(row[col]):
                if not in_rectangle:
                    in_rectangle = True
                    start_col = col
            else:
                if in_rectangle:
                    in_rectangle = False
                    rectangles.append((start_col, col - 1))
        if in_rectangle:
            rectangles.append((start_col, num_cols - 1))
        return rectangles

    def is_closed_rectangle(start_row, end_row, start_col, end_col):
        for col in range(start_col, end_col + 1):
            if not (is_rectangle_char(rows[start_row][col]) and is_rectangle_char(rows[end_row][col])):
                return False
        for row in range(start_row, end_row + 1):
            if not (is_rectangle_char(rows[row][start_col]) and is_rectangle_char(rows[row][end_col])):
                return False
        return True

    rectangles = []
    for row in range(num_rows):
        row_rectangles = find_rectangles_in_row(rows[row])
        for start_col, end_col in row_rectangles:
            for end_row in range(row + 1, num_rows):
                if is_closed_rectangle(row, end_row, start_col, end_col):
                    rectangles.append((row, end_row, start_col, end_col))
                else:
                    break

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
def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    num_cols = len(rows[0]) if num_rows > 0 else 0
    visited = set()
    rectangles = set()

    def mark_rectangle(r, c, dr, dc):
        # Mark the rectangle starting at (r, c) and extending in direction (dr, dc)
        start = (r, c)
        while 0 <= r < num_rows and 0 <= c < num_cols and rows[r][c] in '#█':
            visited.add((r, c))
            r += dr
            c += dc
        return (r - dr, c - dc)

    for r in range(num_rows):
        for c in range(num_cols):
            if (r, c) not in visited and rows[r][c] in '#█':
                # Find the bottom-right corner of the rectangle
                br_r, br_c = mark_rectangle(r, c, 0, 1)  # Horizontal
                br_r, br_c = mark_rectangle(br_r, br_c, 1, 0)  # Vertical
                rectangles.add((r, c, br_r, br_c))

    return len(rectangles)

# The ASCII grid as a string
ascii_grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                  #######       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #     #       
                                                                  #   ##█##     
            ######################################################█###█#█#█#### 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #                                                     #   # # #   # 
            #               #######█████████████████##############█###█#█##   # 
            #         ######█##    #               #  ####        #   # ###   # 
            #         #     # #    #               #  #  #        #   # ###   # 
            ##########█#####█#█####█###############█##█##█########█###█#███#### 
                      #     # #  ##█####           #  #  #        #   # ###     
                      #     # #  # #   #           #  #  #        #   # ###     
                      #     # #  # #   #           #  #  #        #   # ###     
                      #     # #  # #   #           #  #  #        #   # ###     
      ################█#####█#█##█#█###█#          #  #  #        #   # ###     
      #               #     ##█##█#█###██###########  #  #        #   # ###     
      #               #       #  # ####██#############█##█########█###█#█##     
      #               #       #  ########             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        #   # # #     
      #               #       #         #             #  #        ####█## #     
      ################█#######█##########             #  #            #   #     
                      #       #                       #  #            #   #     
                  ####█#######█##########             #  #            #   #     
                  #   #########         #             #  #            #####     
                  #                     #             #  #                      
                  #                     #             #  #                      
                  #                     #             #  #                      
                  #                     #             #  #                      
                  #                     #             #  #                      
                  #                     #             #  #                      
                  #                     #             ####                      
                  #                     #                                       
                  #                     #                                       
                  #                     #                                       
                  #                     #                                       
                  #                     #                                       
                  #                     #                                       
                  #                     #                                       
                  #                     #                                       
                  #                     #                                       
                  #######################                                       
"""

# Count the rectangles
rectangle_count = count_rectangles(ascii_grid)
print(rectangle_count)
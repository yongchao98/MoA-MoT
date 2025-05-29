def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    num_cols = len(rows[0]) if num_rows > 0 else 0
    visited = set()
    rectangles = set()

    def find_rectangle(r, c):
        # Find the bottom-right corner of the rectangle starting at (r, c)
        br_r, br_c = r, c
        while br_c < num_cols and rows[r][br_c] in '#█':
            br_c += 1
        br_c -= 1
        while br_r < num_rows and all(rows[br_r][cc] in '#█' for cc in range(c, br_c + 1)):
            br_r += 1
        br_r -= 1
        return br_r, br_c

    for r in range(num_rows):
        for c in range(num_cols):
            if (r, c) not in visited and rows[r][c] in '#█':
                br_r, br_c = find_rectangle(r, c)
                rectangles.add((r, c, br_r, br_c))
                for rr in range(r, br_r + 1):
                    for cc in range(c, br_c + 1):
                        visited.add((rr, cc))

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
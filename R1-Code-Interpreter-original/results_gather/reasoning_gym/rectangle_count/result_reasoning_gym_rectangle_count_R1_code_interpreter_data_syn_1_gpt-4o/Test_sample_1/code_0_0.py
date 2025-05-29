def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = len(rows[0]) if height > 0 else 0
    visited = set()
    rectangles = set()

    def find_rectangle(x, y):
        if (x, y) in visited or rows[y][x] not in '#█':
            return None
        # Find the bottom-right corner of the rectangle
        max_x, max_y = x, y
        while max_x + 1 < width and rows[y][max_x + 1] in '#█':
            max_x += 1
        while max_y + 1 < height and all(rows[max_y + 1][i] in '#█' for i in range(x, max_x + 1)):
            max_y += 1
        # Mark all parts of the rectangle as visited
        for i in range(x, max_x + 1):
            for j in range(y, max_y + 1):
                visited.add((i, j))
        return (x, y, max_x, max_y)

    for y in range(height):
        for x in range(width):
            if (x, y) not in visited and rows[y][x] in '#█':
                rect = find_rectangle(x, y)
                if rect:
                    rectangles.add(rect)

    return len(rectangles)

# The ASCII grid as a string
ascii_grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                        ############            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #          #            
                                                        #    ######█############
                                                        #    #     #           #
                                                        #    #     #           #
                                                        #    #     #           #
                                                        #    #     #           #
                                                        #    #     #           #
                     ###################################█####█##   #           #
                     #                                  #    # #   #           #
                     #                                  #    # #   #   ########█
                     #                                  #    # #   #   #       █
                     #                                  #    # #   #   #       █
                     #                                  #    # #   #   #       █
      ###############█##################################█####█##   #   #       █
      ###############█############################      #    ###   #   #       █
      ##             #                           #      #    ###   #   #       █
      ##             #                           #      #    ###   #   #       █
      ##             #                           #      #    ###   #   #       █
      ##             #                           #      #    ###   #   ########█
      ##             #                           #      #    ###   #           #
      ##             #                           #      #    ###   #           #
      ##             #                           #      #    ###   #           #
      ##             #                           #      #    ###   #           #
      ##             #                  ###      #      #    ###   #           #
      ##             #                  # #      #      #    ###   #           #
      ##             #                  # #      #      #    ###   #           #
      ##             ###################█#█######█######█####██#   #           #
      ##                                # #      #      #####██#####           #
      ##                                # #      #           ##                #
      ##################################█#█#######           ##                #
      #                                 # #                  ##                #
      #                                 # #                  ##                #
      #                                 # #                  ##                #
      #                                 # #                  ##                #
      #                                 # #                  ##                #
      ##################################█#█##################█#                #
                                        # #                  #                 #
                                        # #                  #                 #
                                        # #                  #                 #
                                        # #                  #                 #
                                        # #                  #                 #
                                        # #                  #                 #
                                        # #                  #                 #
                                        ##█################  #                 #
                                        #█#         ######█##███################
                                         #          #     #    #########        
                                         #          #     #    ##      #        
                                         #          #     #    ##      #        
                                         #          #     #    ##      #        
                                         ###########█######    ##      #        
                                                    #          ##      #        
                                                    #          ##      #        
                                                    #          ##      #        
                                                    #############      #        
                                                                ########        
                                                                                
                                                                                
                                                                                
"""

# Count the rectangles
rectangle_count = count_rectangles(ascii_grid)
print(rectangle_count)
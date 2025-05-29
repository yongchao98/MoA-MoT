def count_rectangles(grid):
    lines = grid.strip().split('\n')
    height = len(lines)
    width = len(lines[0])
    
    # First, let's count distinct rectangles
    distinct_rectangles = 0
    
    # Find top-left corners and validate rectangles
    for y in range(height):
        for x in range(width):
            # Check if this could be a top-left corner
            if lines[y][x] == '#':
                # Verify it's actually a top-left corner
                if (y == 0 or lines[y-1][x] not in '#█') and (x == 0 or lines[y][x-1] not in '#█'):
                    distinct_rectangles += 1

    # Count overlapping points (█)
    overlapping_points = sum(line.count('█') for line in lines)
    
    # Each █ represents an additional rectangle
    total_rectangles = distinct_rectangles + overlapping_points
    
    print(total_rectangles)

# Input grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                ###########     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                             ###################################█##       #     
                             #                                  # #       #     
                             #                                  # #       #     
                             #                                  # #       #     
                             #                                  # #       #     
                             #                                  ##█########     
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #        #####
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             ######################################        #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                      #####█###█
                                                                      #    #   █
                                                                      #    ####█
                                                                      #        #
                                                                      #        #
                                                                      #        #
                                                                      #        #
                                                                      ##########
                                                                                
                                                                                """

count_rectangles(grid)
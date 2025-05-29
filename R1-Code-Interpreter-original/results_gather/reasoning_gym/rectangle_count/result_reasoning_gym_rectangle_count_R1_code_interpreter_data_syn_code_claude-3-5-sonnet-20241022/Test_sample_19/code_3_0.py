def count_rectangles(grid):
    lines = grid.strip().split('\n')
    height = len(lines)
    width = len(lines[0])
    
    def find_rectangle_corners():
        corners = []
        # Find all top-left corners
        for y in range(height):
            for x in range(width):
                if lines[y][x] in '#█':
                    # Check if it's a top-left corner
                    is_top = y == 0 or lines[y-1][x] not in '#█'
                    is_left = x == 0 or lines[y][x-1] not in '#█'
                    if is_top and is_left:
                        # Find bottom-right corner
                        end_x = x
                        while end_x < width and lines[y][end_x] in '#█':
                            end_x += 1
                        end_x -= 1
                        
                        end_y = y
                        while end_y < height and lines[end_y][x] in '#█':
                            end_y += 1
                        end_y -= 1
                        
                        # Verify rectangle
                        valid = True
                        for cy in range(y, end_y + 1):
                            for cx in range(x, end_x + 1):
                                if cy in (y, end_y) or cx in (x, end_x):
                                    if lines[cy][cx] not in '#█':
                                        valid = False
                                        break
                        if valid:
                            corners.append((y, x, end_y, end_x))
    
        return corners

    def count_overlaps():
        count = 0
        for y in range(height):
            for x in range(width):
                if lines[y][x] == '█':
                    count += 1
        return count

    base_rectangles = len(find_rectangle_corners())
    overlaps = count_overlaps()
    
    # Each '█' represents an additional rectangle
    total = base_rectangles + overlaps
    
    print(total)

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
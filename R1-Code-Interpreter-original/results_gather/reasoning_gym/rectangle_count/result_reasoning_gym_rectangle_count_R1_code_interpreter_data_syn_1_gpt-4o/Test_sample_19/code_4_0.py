def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    num_cols = len(rows[0]) if num_rows > 0 else 0
    visited = [[False] * num_cols for _ in range(num_rows)]
    rectangle_count = 0

    def mark_visited(r, c, height, width):
        for i in range(r, r + height):
            for j in range(c, c + width):
                visited[i][j] = True

    for r in range(num_rows):
        for c in range(num_cols):
            if (rows[r][c] == '#' or rows[r][c] == '█') and not visited[r][c]:
                # Determine the dimensions of the rectangle
                width = 0
                while c + width < num_cols and (rows[r][c + width] == '#' or rows[r][c + width] == '█'):
                    width += 1
                height = 0
                while r + height < num_rows and (rows[r + height][c] == '#' or rows[r + height][c] == '█'):
                    height += 1

                # Mark the rectangle as visited
                mark_visited(r, c, height, width)

                # Count the rectangle
                rectangle_count += 1

    return rectangle_count

# The given ASCII grid
ascii_grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
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

# Count the rectangles
rectangle_count = count_rectangles(ascii_grid)
print(rectangle_count)
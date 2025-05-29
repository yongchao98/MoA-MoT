def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = len(rows[0]) if height > 0 else 0

    def is_rectangle_char(c):
        return c == '#' or c == '█'

    def find_rectangle_end(r, c):
        # Find the bottom-right corner of the rectangle starting at (r, c)
        if not is_rectangle_char(rows[r][c]):
            return None

        # Find the width of the rectangle
        rect_width = 0
        while c + rect_width < width and is_rectangle_char(rows[r][c + rect_width]):
            rect_width += 1

        # Find the height of the rectangle
        rect_height = 0
        while r + rect_height < height and all(is_rectangle_char(rows[r + rect_height][c + i]) for i in range(rect_width)):
            rect_height += 1

        return (r + rect_height - 1, c + rect_width - 1)

    visited = set()
    rectangle_count = 0

    for r in range(height):
        for c in range(width):
            if (r, c) not in visited and is_rectangle_char(rows[r][c]):
                end = find_rectangle_end(r, c)
                if end:
                    rectangle_count += 1
                    # Mark all cells in this rectangle as visited
                    for i in range(r, end[0] + 1):
                        for j in range(c, end[1] + 1):
                            visited.add((i, j))

    return rectangle_count

# Define the grid
grid = """
                                                                                
                                                                                
                                                                                
                                               #######                          
                                               #     #                          
                                               #     #                          
                                               #     #                          
                                               #     #                          
                                               #     #                          
                                               #     #                          
                                               #     #                          
                                               #     #                          
                                               #     #             #############
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #           #
                                               #     #             #      ######
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                                               #     #             #      #   ##
                ###############################█#####█####         #      #   ##
                #                              #     #   #         #      #   ##
                #                              #     #   #         #      #   ##
                #                              #     #   #         #######█###█#
                #                              #     #   #                #   # 
                #                              #     #   #                #   # 
                #                              #     #   #                #   # 
                #                              #     #   ####             #   # 
                #                              #     #   ## #             #   # 
                ###############################█#####█##### #             #   # 
                                               #     #    # #             #   # 
                                               #     #    # #             #   # 
                                               #     #    # #             #   # 
                                               #     #    # #             ##### 
                                               #     #    # #                   
                                             ##█#    #    # #                   
                                             # ##    #    # #                   
                                             ##█#    #    # #                   
                                               #     #####█#█####               
                                               #     ##   # #   #               
                                               #     ##   # #   #               
                                               #     ##   # #   #               
                                               #     ##   # #   #               
                                               #     ##   # #   #               
                                               #     #█###███###█##             
                                               #     #█         # #             
                                               #     #█         # #             
                                               #     #█         # #   #####     
                                               #     #█         # #   #   #     
                                               #     #█         # #   #   #     
                                               #######█         # #   #   #     
                                                      █         # #   #   #     
                                                      █         # #   #   #     
                                                      █         # #   #   #     
                                                      █         # #   #   #     
                                                      █         # #   #   #     
                                                      █         # #   #   #     
                                                      █         # #   #   #     
                                                      █#########█##   #   #     
                                                      #         #     #   #     
                                                      ###########     #####     
"""

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
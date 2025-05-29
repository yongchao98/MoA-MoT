def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    # Find all '#' and '█' positions
    rectangles = []
    height = len(grid)
    width = len(grid[0])
    
    # Function to check if position is part of existing rectangle
    def is_part_of_existing(x, y):
        for rect in rectangles:
            if (rect[0] <= x <= rect[2] and rect[1] <= y <= rect[3]):
                return True
        return False
    
    # Find rectangles
    for y in range(height):
        for x in range(width):
            if (grid[y][x] in '#█') and not is_part_of_existing(x, y):
                # Found top-left corner of new rectangle
                start_x, start_y = x, y
                
                # Find width
                end_x = start_x
                while end_x < width and grid[y][end_x] in '#█':
                    end_x += 1
                end_x -= 1
                
                # Find height
                end_y = start_y
                while end_y < height and grid[end_y][x] in '#█':
                    end_y += 1
                end_y -= 1
                
                # Verify rectangle
                is_rectangle = True
                for cy in range(start_y, end_y + 1):
                    for cx in range(start_x, end_x + 1):
                        if cx in (start_x, end_x) or cy in (start_y, end_y):
                            if grid[cy][cx] not in '#█':
                                is_rectangle = False
                                break
                
                if is_rectangle:
                    rectangles.append((start_x, start_y, end_x, end_y))
    
    # Count overlapping rectangles
    overlap_count = 0
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                count = sum(1 for rect in rectangles if rect[0] <= x <= rect[2] and rect[1] <= y <= rect[3])
                if count == 2:  # If position is part of two rectangles
                    overlap_count += 1
                    break  # We only need to find one overlapping point to confirm

    # Total rectangles = number of distinct rectangles + number of additional rectangles from overlaps
    total_rectangles = len(rectangles)
    if overlap_count > 0:
        total_rectangles += 1  # Add one more rectangle for overlapping area

    print(total_rectangles)

# Your ASCII grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                         ####################################################   
                         #                                                  #   
                         #                                                  #   
                         #                                                  #   
                         #                                                  #   
                         #                                            ##### #   
                         #                                            #   # #   
                         #                                            #   # #   
                         #                                            #   # #   
                         #                                            ##### #   
                         #                                                  #   
                         #                                                  #   
                         #                                                  #   
                         #                                                  #   
                         #                                                  #   
                         #                                                  #   
                         ####################################################   
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                          ######################################
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                    ######█#################                   #
                                    █#####█#######         #                   #
                                    █     #      #         #                   #
                                    █#####█######█##########                   #
                                    #     #      #                             #
                                    #     #      #                             #
                                    #     #      #                             #
                                    #     #      #                             #
                                    #     #      #                             #
                                    #     #      #                             #
                                    ######█#######                             #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          #                                    #
                                          ######################################"""

find_rectangles(grid)
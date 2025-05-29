def count_rectangles(grid):
    # Convert the grid into a list of lines
    lines = grid.strip().split('\n')
    
    # Initialize variables to store the number of rectangles
    single_rectangles = 0
    overlapping_rectangles = 0
    
    # Scan the grid for rectangles
    for i in range(len(lines)):
        for j in range(len(lines[i])):
            if lines[i][j] == '#':
                # Check if this is the top-left corner of a rectangle
                if (i == 0 or lines[i-1][j] not in ['#', '█']) and (j == 0 or lines[i][j-1] not in ['#', '█']):
                    single_rectangles += 1
            elif lines[i][j] == '█':
                # Check if this is the top-left corner of an overlapping rectangle
                if (i == 0 or lines[i-1][j] != '█') and (j == 0 or lines[i][j-1] != '█'):
                    overlapping_rectangles += 1
    
    # Total rectangles is the sum of single and overlapping rectangles
    total_rectangles = single_rectangles + overlapping_rectangles
    return total_rectangles

# Define the grid
grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                      ###       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                            ##########################█#█       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         # █       
                                            #                         ##█       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #                           #       
                                            #############################       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
"""

# Count the rectangles
total_rectangles = count_rectangles(grid)
print(total_rectangles)
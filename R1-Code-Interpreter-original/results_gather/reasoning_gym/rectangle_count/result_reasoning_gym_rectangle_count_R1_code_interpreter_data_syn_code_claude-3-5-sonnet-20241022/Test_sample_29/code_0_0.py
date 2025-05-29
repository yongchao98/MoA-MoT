def is_corner(grid, i, j):
    return grid[i][j] in ['#', '█']

def verify_rectangle(grid, top, left, bottom, right):
    # Check if all edges are made of # or █
    # Top and bottom edges
    for j in range(left, right + 1):
        if not grid[top][j] in ['#', '█'] or not grid[bottom][j] in ['#', '█']:
            return False
    
    # Left and right edges
    for i in range(top, bottom + 1):
        if not grid[i][left] in ['#', '█'] or not grid[i][right] in ['#', '█']:
            return False
    
    return True

def count_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Convert input to 2D array
    grid = [list(row) for row in grid.split('\n')]
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if not is_corner(grid, top, left):
                continue
                
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if (is_corner(grid, top, left) and 
                        is_corner(grid, top, right) and 
                        is_corner(grid, bottom, left) and 
                        is_corner(grid, bottom, right)):
                        
                        if verify_rectangle(grid, top, left, bottom, right):
                            rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# Input grid
grid = """                                                                                
                                                                                
                                                                                
                                     #######################                    
                                     #                     #                    
                                     #                     #                    
                                     #                     #                    
                ##############       #                     #                    
                #            #       #       ##############█##########          
                #            #       #       #             #         #          
                #            #       #       #             #         #          
################█############█#######█#######█#############█#######  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       #             #      #  #          
#               #            #       #       # ############█######█##█####      
#               #            #       #       # #           #      #  #   #      
#               #            #       #       # #           #      #  #   #      
#               #            #       #       # #           #      #  #   #      
#               #            #       #       # #           #      #  #   #      
#               #      ######█#######█#######█#█###########█#     #  #   #      
#               #      #     #       #       # #           ##     #  #   #      
#               #      #     #       #       # #           ##     #  #   #      
#               #      #     #       #       # #           ##     #  #   #      
#               #      #     #       #       # #           ##     #  #   #      
#               #      #     #       #       # #           ##     #  #   #      
#               #      #     #       #       # #           ##     #  #   #      
#               #      #     #       #       # #           ##     #  #   #      
################█######█#####█#######█#######█#█###########██######  #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      #     #       #       # #           ##        #   #      
                #      ######█#######█#######█#█###########█#        #   #      
                #            #       #       # #           #         #   #      
                █############█#######█#######█#█###########█#####    #   #      
                █            #       #       # #           #    #    #   #      
                █            #       #       # #           #    #    #   #      
                █            #       #       # #           #    #    #   #      
                █            #       #     ##█##           #    #    #   #      
                █            #       #     # ###           #    #    #   #      
                █            #       #     # ###           #    #    #   #      
                █            #       #     # ###           #    #    #   #      
                █            #       #     # ###           #    #    #   #      
                █            #       #     # ###           #    #    #   #      
                █            #       #     ##███###########█####█#####   #      
                █############█#######█#########█###########█#####        #      
                #            #       #         #           #             #      
                #            #       #         #           #             #      
                #            #       #         #           #             #      
                #            #       #         #           #             #      
                #            #       ##########█############             #      
                #            #                 #                         #      
                #            #                 #                         #      
                #            #                 #                         #      
                #            #                 #                         #      
                #            #                 #                         #      
                #            #                 #                         #      
                #            #                 ###########################      
                #            #                                                  
                #            #                                                  
                #            #                                                  
                ##############                                                  
                                                                                """

print(count_rectangles(grid))
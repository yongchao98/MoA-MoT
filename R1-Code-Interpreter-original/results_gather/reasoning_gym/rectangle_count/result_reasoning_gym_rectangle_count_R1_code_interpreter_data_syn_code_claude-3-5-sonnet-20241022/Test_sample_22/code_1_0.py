def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    # Find dimensions
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    # Find all border points
    border_points = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                border_points.append((y, x))
    
    # Function to check if a rectangle exists between two diagonal points
    def is_rectangle(y1, x1, y2, x2):
        if y1 >= y2 or x1 >= x2:
            return False
            
        # Check all four edges
        for x in range(x1, x2 + 1):
            if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                return False
        for y in range(y1, y2 + 1):
            if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                return False
        return True
    
    # Count rectangles
    rectangles = set()
    for i, (y1, x1) in enumerate(border_points):
        for y2, x2 in border_points[i+1:]:
            if y2 > y1 and x2 > x1:  # Only check points that could form a rectangle
                if is_rectangle(y1, x1, y2, x2):
                    rectangles.add((y1, x1, y2, x2))
    
    print(len(rectangles))

# Create the grid string (your provided grid here)
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                          ###########################           
                                          #                         #           
                                          #                   ######█########   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #             ######█#####█########   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     ######█######█#   
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             ############█#######    
                                          #                         #           
                                          #                         #           
                                          #                         #           
                                          ###########################           
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)
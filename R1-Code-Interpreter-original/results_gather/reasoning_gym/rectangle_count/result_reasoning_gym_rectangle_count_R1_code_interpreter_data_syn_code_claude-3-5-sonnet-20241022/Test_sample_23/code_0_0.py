def find_rectangles(grid):
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        
        # Check left and right edges
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        
        return True
    
    rectangles = 0
    
    # Find all possible rectangles by checking each possible corner pair
    for top in range(rows):
        for left in range(cols):
            if grid[top][left] not in '#█':
                continue
                
            for bottom in range(top + 1, rows):
                for right in range(left + 1, cols):
                    if grid[bottom][right] not in '#█':
                        continue
                        
                    if is_valid_rectangle(top, left, bottom, right):
                        rectangles += 1
    
    return rectangles

# Create the grid from the ASCII art
grid = []
empty_line_count = 0
non_empty_lines = []

# Read all non-empty lines and store them
for line in """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                       ######################                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       # ############       #                                   
                       # #          #       #                                   
                       # #          #       #                                   
                       # #          #       #                                   
                       # #          #       #                                   
                       # #          #       #                                   
                       # #          #       #                                   
                       # ############       #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       #                    #                                   
                       ######################                                   
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
""".strip('\n').split('\n'):
    if line.strip():
        non_empty_lines.append(line)

grid = non_empty_lines
result = find_rectangles(grid)
print(result)
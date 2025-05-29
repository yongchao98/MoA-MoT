def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line.rstrip() for line in grid.split('\n')]
    
    # Find the actual grid boundaries by removing empty lines
    start_y = 0
    end_y = len(lines)
    while start_y < end_y and not any(c in '#█' for c in lines[start_y]):
        start_y += 1
    while end_y > start_y and not any(c in '#█' for c in lines[end_y-1]):
        end_y -= 1
    
    if start_y >= end_y:
        return 0
    
    # Get actual grid
    grid = lines[start_y:end_y]
    height = len(grid)
    width = max(len(row) for row in grid)
    
    # Standardize grid width
    grid = [row.ljust(width) for row in grid]
    
    count = 0
    visited = set()
    
    def trace_rectangle(y, x):
        if (y, x) in visited:
            return False
        
        # Find right boundary
        right = x
        while right < width and grid[y][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom boundary
        bottom = y
        while bottom < height and grid[bottom][x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        is_rectangle = True
        for cy in range(y, bottom + 1):
            if grid[cy][x] not in '#█' or grid[cy][right] not in '#█':
                is_rectangle = False
                break
        for cx in range(x, right + 1):
            if grid[y][cx] not in '#█' or grid[bottom][cx] not in '#█':
                is_rectangle = False
                break
                
        if is_rectangle:
            # Mark all border points as visited
            for cy in range(y, bottom + 1):
                visited.add((cy, x))
                visited.add((cy, right))
            for cx in range(x, right + 1):
                visited.add((y, cx))
                visited.add((bottom, cx))
            return True
        return False
    
    # Scan for rectangles
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                if trace_rectangle(y, x):
                    count += 1
    
    print(count)

# Your grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                       #########################################
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                       #
                                       #                                    ####
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #                                    # ##
                                       #####################################█#█#
                                                                            # # 
                                                                            # # 
                                                                            # # 
                                                                            # # 
                                                                            # # 
                                                                            # # 
                                                                            ### 
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)
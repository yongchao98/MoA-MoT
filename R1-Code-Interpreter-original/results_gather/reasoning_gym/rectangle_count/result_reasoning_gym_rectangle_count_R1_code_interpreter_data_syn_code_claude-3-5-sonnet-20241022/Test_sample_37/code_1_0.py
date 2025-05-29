def is_complete_rectangle(grid, start_row, start_col):
    # Find width of top line
    width = 0
    for col in range(start_col, len(grid[start_row])):
        if grid[start_row][col] == '#':
            width += 1
        else:
            break
    
    if width < 2:  # Need at least 2 characters for a rectangle
        return False
    
    # Find height by checking first column
    height = 0
    for row in range(start_row, len(grid)):
        if row < len(grid) and grid[row][start_col] == '#':
            height += 1
        else:
            break
    
    if height < 2:  # Need at least 2 characters for a rectangle
        return False
    
    # Check if it forms a complete rectangle
    # Check top and bottom rows
    for col in range(start_col, start_col + width):
        if (grid[start_row][col] != '#' or 
            grid[start_row + height - 1][col] != '#'):
            return False
    
    # Check left and right columns
    for row in range(start_row, start_row + height):
        if (grid[row][start_col] != '#' or 
            grid[row][start_col + width - 1] != '#'):
            return False
    
    # Check if this is the top-left corner
    if (start_row > 0 and grid[start_row-1][start_col] == '#'):
        return False
    if (start_col > 0 and grid[start_row][start_col-1] == '#'):
        return False
    
    return True

# Initialize the grid
grid_str = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
            ########################                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            ########################                                            
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

# Convert string to list of lines
grid = grid_str.splitlines()

# Count rectangles
count = 0
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == '#' and is_complete_rectangle(grid, i, j):
            count += 1

print(count)
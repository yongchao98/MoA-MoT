def count_rectangles(grid):
    # Convert the grid into a list of lines
    lines = grid.strip().split('\n')
    
    # Initialize variables
    rectangle_count = 0
    in_rectangle = False
    
    # Iterate over each line
    for line in lines:
        # Check if the line contains a rectangle boundary
        if '###' in line:
            if not in_rectangle:
                rectangle_count += 1
                in_rectangle = True
        else:
            in_rectangle = False
    
    return rectangle_count

# Define the grid
grid = """
                                                                                
                                                                                
                                                                                
                                                                          ###   
                                                                          # #   
                                                                          # #   
                                                                          # #   
                                                                          # #   
                                                                          ###   
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
"""

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
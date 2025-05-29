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
                # Start of a new rectangle
                in_rectangle = True
        elif in_rectangle and line.strip() == '':
            # End of a rectangle
            rectangle_count += 1
            in_rectangle = False
    
    # If the last line was part of a rectangle, count it
    if in_rectangle:
        rectangle_count += 1
    
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
def find_rectangles():
    rectangles = 0
    grid = []
    
    # Read only non-empty lines
    try:
        while True:
            line = input()
            if '#' in line or '█' in line:
                grid.append(line)
    except EOFError:
        pass
    
    if not grid:
        print("<<<0>>>")
        return
    
    # Find first and last row with content
    start_row = 0
    end_row = len(grid) - 1
    
    # Process each rectangle by finding its boundaries
    for row in range(start_row, end_row + 1):
        line = grid[row]
        start_col = 0
        while start_col < len(line):
            # Find start of potential rectangle
            while start_col < len(line) and line[start_col] not in '#█':
                start_col += 1
            if start_col >= len(line):
                break
                
            # Find end of potential rectangle
            end_col = start_col
            while end_col < len(line) and line[end_col] in '#█':
                end_col += 1
            
            # If we found a horizontal line, check if it's part of a rectangle
            if end_col - start_col > 1:
                # Look for matching bottom line
                for bottom_row in range(row + 1, end_row + 1):
                    bottom_line = grid[bottom_row]
                    
                    # Quick check if this could be bottom of rectangle
                    if all(bottom_line[x] in '#█' for x in [start_col, end_col-1]):
                        # Verify vertical lines
                        is_rectangle = True
                        for check_row in range(row + 1, bottom_row):
                            if (grid[check_row][start_col] not in '#█' or 
                                grid[check_row][end_col-1] not in '#█'):
                                is_rectangle = False
                                break
                        
                        if is_rectangle and all(bottom_line[x] in '#█' 
                                              for x in range(start_col, end_col)):
                            rectangles += 1
            
            start_col = end_col + 1
    
    print(f"<<<{rectangles}>>>")

# Run the analysis
find_rectangles()
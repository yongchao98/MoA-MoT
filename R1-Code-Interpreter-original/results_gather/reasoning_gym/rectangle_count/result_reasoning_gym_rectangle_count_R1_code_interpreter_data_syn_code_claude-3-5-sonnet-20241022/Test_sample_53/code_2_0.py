def scan_grid():
    # Initialize empty grid with fixed size
    grid = [[' ' for _ in range(100)] for _ in range(100)]
    row = 0
    
    # Read input more efficiently
    try:
        while True:
            line = input()
            for col, char in enumerate(line):
                grid[row][col] = char
            row += 1
    except EOFError:
        pass
    
    rectangles = 0
    height = row
    width = 100  # Maximum width
    
    # Store horizontal lines (y, x1, x2)
    h_lines = []
    
    # Find horizontal lines
    for y in range(height):
        start = None
        for x in range(width):
            if grid[y][x] in '#█':
                if start is None:
                    start = x
            elif start is not None:
                if x - start > 1:  # Minimum width of 2
                    h_lines.append((y, start, x-1))
                start = None
        if start is not None and width - start > 1:
            h_lines.append((y, start, width-1))
    
    # Check each pair of horizontal lines
    for i, (y1, x1_start, x1_end) in enumerate(h_lines):
        for y2, x2_start, x2_end in h_lines[i+1:]:
            if y2 <= y1:
                continue
            
            # Find overlapping x-range
            x_start = max(x1_start, x2_start)
            x_end = min(x1_end, x2_end)
            
            if x_end - x_start < 1:
                continue
            
            # Check vertical lines
            valid = True
            for x in (x_start, x_end):
                for y in range(y1+1, y2):
                    if grid[y][x] not in '#█':
                        valid = False
                        break
                if not valid:
                    break
            
            if valid:
                rectangles += 1
    
    print(f"<<<{rectangles}>>>")

# Run the scanner
scan_grid()
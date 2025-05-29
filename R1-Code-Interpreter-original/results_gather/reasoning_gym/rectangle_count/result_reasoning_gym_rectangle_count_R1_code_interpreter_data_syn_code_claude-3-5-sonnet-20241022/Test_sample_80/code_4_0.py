def quick_count_rectangles(grid):
    # Process only non-empty lines to reduce data
    non_empty_lines = []
    start_line = 0
    for i, line in enumerate(grid):
        if any(c in '#█' for c in line):
            if not non_empty_lines:
                start_line = i
            non_empty_lines.append(line)
    
    if not non_empty_lines:
        return 0
    
    # Count obvious rectangles
    count = 0
    width = len(non_empty_lines[0])
    height = len(non_empty_lines)
    
    # Simple rectangle detection
    for i in range(height):
        for j in range(width):
            if non_empty_lines[i][j] in '#█':
                # Check if start of rectangle
                if (j == 0 or non_empty_lines[i][j-1] not in '#█') and (i == 0 or non_empty_lines[i-1][j] not in '#█'):
                    # Find right edge
                    right = j
                    while right < width and non_empty_lines[i][right] in '#█':
                        right += 1
                    right -= 1
                    
                    # Find bottom edge
                    bottom = i
                    while bottom < height and non_empty_lines[bottom][j] in '#█':
                        bottom += 1
                    bottom -= 1
                    
                    # Quick validate
                    if bottom > i and right > j:
                        is_valid = True
                        has_overlap = False
                        
                        # Check only corners for validation
                        corners = [
                            (i, j), (i, right),
                            (bottom, j), (bottom, right)
                        ]
                        
                        for y, x in corners:
                            if non_empty_lines[y][x] not in '#█':
                                is_valid = False
                                break
                            if non_empty_lines[y][x] == '█':
                                has_overlap = True
                        
                        if is_valid:
                            count += 1
                            if has_overlap:
                                count += 1
    
    return count

# Read only non-empty lines
grid = []
empty_line_count = 0
for _ in range(80):  # Fixed number of iterations
    try:
        line = input()
        if any(c in '#█' for c in line):
            grid.append(line)
            empty_line_count = 0
        else:
            empty_line_count += 1
            if empty_line_count > 5:  # Stop if too many empty lines
                break
    except EOFError:
        break

print(f"<<<{quick_count_rectangles(grid)}>>>")
def find_rectangles():
    # Store horizontal and vertical line segments as tuples
    h_lines = []  # (y, x1, x2)
    v_lines = []  # (x, y1, y2)
    
    # Read and process input efficiently
    y = 0
    try:
        while True:
            line = input().rstrip()
            if not line:
                y += 1
                continue
                
            # Find horizontal lines
            start = None
            for x, char in enumerate(line):
                if char in '#█':
                    if start is None:
                        start = x
                elif start is not None:
                    if x - start > 1:  # Minimum width check
                        h_lines.append((y, start, x - 1))
                    start = None
            if start is not None and len(line) - start > 1:
                h_lines.append((y, start, len(line) - 1))
            y += 1
    except EOFError:
        pass
    
    if not h_lines:
        return 0
    
    # Convert horizontal lines to vertical lines
    max_width = max(x2 for _, x1, x2 in h_lines) + 1
    vertical_points = [[] for _ in range(max_width)]
    
    # Map points to columns
    for y, x1, x2 in h_lines:
        for x in range(x1, x2 + 1):
            vertical_points[x].append(y)
    
    # Find vertical lines
    for x, points in enumerate(vertical_points):
        if not points:
            continue
        points.sort()
        start = None
        for y in points:
            if start is None:
                start = y
            elif y - points[points.index(y) - 1] > 1:
                if y - start > 1:  # Minimum height check
                    v_lines.append((x, start, points[points.index(y) - 1]))
                start = y
        if start is not None and points[-1] - start > 1:
            v_lines.append((x, start, points[-1]))
    
    # Count rectangles using line intersections
    count = 0
    for i, (y1, left1, right1) in enumerate(h_lines):
        for y2, left2, right2 in h_lines[i + 1:]:
            if y2 <= y1:
                continue
            
            # Find overlapping x-range
            left = max(left1, left2)
            right = min(right1, right2)
            
            if left >= right:
                continue
            
            # Check for vertical lines at boundaries
            has_left = False
            has_right = False
            
            for x, top, bottom in v_lines:
                if x == left and top <= y1 and bottom >= y2:
                    has_left = True
                elif x == right and top <= y1 and bottom >= y2:
                    has_right = True
                if has_left and has_right:
                    count += 1
                    break
    
    return count

print(find_rectangles())
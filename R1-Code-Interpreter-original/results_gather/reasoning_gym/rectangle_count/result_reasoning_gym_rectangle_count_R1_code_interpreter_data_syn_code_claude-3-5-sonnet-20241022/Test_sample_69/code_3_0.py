def process_grid():
    # Process input and create compressed grid representation
    points = set()
    min_x, max_x, min_y, max_y = float('inf'), 0, float('inf'), 0
    y = 0
    
    try:
        while True:
            line = input()
            if not any(c in '#█' for c in line):
                y += 1
                continue
            
            for x, c in enumerate(line):
                if c in '#█':
                    points.add((x, y))
                    min_x = min(min_x, x)
                    max_x = max(max_x, x)
                    min_y = min(min_y, y)
                    max_y = max(max_y, y)
            y += 1
    except EOFError:
        pass
    
    return points, min_x, max_x, min_y, max_y

def count_rectangles():
    points, min_x, max_x, min_y, max_y = process_grid()
    if not points:
        return 0
    
    # Get sorted x and y coordinates
    x_coords = sorted({x for x, _ in points})
    y_coords = sorted({y for _, y in points})
    
    count = 0
    # Check each possible top-left corner
    for i, top in enumerate(y_coords[:-1]):
        for j, left in enumerate(x_coords[:-1]):
            if (left, top) not in points:
                continue
            
            # For each possible bottom-right corner
            for bottom in y_coords[i+1:]:
                for right in x_coords[j+1:]:
                    # Quick check for corners
                    if not all((left, top) in points,
                             (right, top) in points,
                             (left, bottom) in points,
                             (right, bottom) in points):
                        continue
                    
                    # Check top and bottom edges
                    valid = True
                    for x in x_coords:
                        if left < x < right:
                            if (x, top) not in points or (x, bottom) not in points:
                                valid = False
                                break
                    
                    if not valid:
                        continue
                    
                    # Check left and right edges
                    for y in y_coords:
                        if top < y < bottom:
                            if (left, y) not in points or (right, y) not in points:
                                valid = False
                                break
                    
                    if valid:
                        count += 1
    
    return count

print(count_rectangles())
def find_rectangles(grid):
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    rectangles = set()
    
    def is_valid_char(c):
        return c in ['#', '█']
    
    def trace_rectangle(start_r, start_c):
        # Find the width of the rectangle
        width = 0
        for c in range(start_c, cols):
            if not is_valid_char(grid[start_r][c]):
                break
            width += 1
            
        if width < 2:  # Too narrow to be a rectangle
            return None
            
        # Find the height of the rectangle
        height = 0
        for r in range(start_r, rows):
            if not is_valid_char(grid[r][start_c]):
                break
            height += 1
            
        if height < 2:  # Too short to be a rectangle
            return None
            
        # Verify the rectangle
        for r in range(start_r, start_r + height):
            for c in range(start_c, start_c + width):
                if not is_valid_char(grid[r][c]):
                    return None
                if (r > start_r and r < start_r + height - 1 and 
                    c > start_c and c < start_c + width - 1):
                    if is_valid_char(grid[r][c]):
                        return None
                    
        return (start_r, start_c, height, width)
    
    # Find all potential rectangles
    for r in range(rows):
        for c in range(cols):
            if is_valid_char(grid[r][c]):
                rect = trace_rectangle(r, c)
                if rect:
                    rectangles.add(rect)
    
    # Count overlapping rectangles
    overlaps = set()
    for r1 in rectangles:
        for r2 in rectangles:
            if r1 != r2:
                r1_points = set((i, j) for i in range(r1[0], r1[0] + r1[2])
                              for j in range(r1[1], r1[1] + r1[3]))
                r2_points = set((i, j) for i in range(r2[0], r2[0] + r2[2])
                              for j in range(r2[1], r2[1] + r2[3]))
                if r1_points & r2_points:  # If there's intersection
                    overlap_points = r1_points & r2_points
                    for point in overlap_points:
                        if grid[point[0]][point[1]] == '█':
                            overlaps.add(tuple(sorted([r1, r2])))
                            
    # Each overlap represents two rectangles that were counted as one
    total_rectangles = len(rectangles) + len(overlaps)
    print(total_rectangles)

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

find_rectangles(grid)
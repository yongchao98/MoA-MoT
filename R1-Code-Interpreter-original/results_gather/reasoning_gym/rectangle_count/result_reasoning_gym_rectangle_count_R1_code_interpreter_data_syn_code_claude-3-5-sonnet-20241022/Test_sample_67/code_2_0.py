def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    def is_valid_char(char):
        return char in '#█'
    
    def validate_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for j in range(left, right + 1):
            if not is_valid_char(lines[top][j]) or not is_valid_char(lines[bottom][j]):
                return False
        
        # Check left and right edges
        for i in range(top, bottom + 1):
            if not is_valid_char(lines[i][left]) or not is_valid_char(lines[i][right]):
                return False
        
        # Check if interior points are empty or █ (for overlaps)
        for i in range(top + 1, bottom):
            for j in range(left + 1, right):
                if lines[i][j] not in ' #█':
                    return False
        
        return True
    
    def find_next_start(start_i, start_j):
        for i in range(start_i, height):
            j_start = start_j if i == start_i else 0
            for j in range(j_start, width):
                if is_valid_char(lines[i][j]):
                    return (i, j)
        return None
    
    rectangles = []
    start_pos = find_next_start(0, 0)
    
    while start_pos:
        top, left = start_pos
        # Find possible right edges
        for right in range(left + 1, width):
            if not is_valid_char(lines[top][right]):
                break
            # Find possible bottom edges
            for bottom in range(top + 1, height):
                if not is_valid_char(lines[bottom][left]):
                    break
                if validate_rectangle(top, left, bottom, right):
                    rectangles.append((top, left, bottom, right))
        
        # Find next starting position
        start_pos = find_next_start(top, left + 1)
    
    # Filter out duplicate rectangles and contained rectangles
    unique_rectangles = set()
    for r in rectangles:
        is_unique = True
        top1, left1, bottom1, right1 = r
        
        for other in rectangles:
            if r == other:
                continue
            top2, left2, bottom2, right2 = other
            
            # Check if this rectangle is contained within another
            if (top2 <= top1 and left2 <= left1 and bottom2 >= bottom1 and right2 >= right1):
                if all(lines[i][j] == '#' for i in [top1, bottom1] for j in range(left1, right1 + 1)) and \
                   all(lines[i][j] == '#' for j in [left1, right1] for i in range(top1, bottom1 + 1)):
                    is_unique = False
                    break
        
        if is_unique:
            unique_rectangles.add(r)
    
    return len(unique_rectangles)

# Your ASCII grid
test_grid = """                                                                                
[... rest of your grid ...]
                                                                                """

result = find_rectangles(test_grid)
print(result)
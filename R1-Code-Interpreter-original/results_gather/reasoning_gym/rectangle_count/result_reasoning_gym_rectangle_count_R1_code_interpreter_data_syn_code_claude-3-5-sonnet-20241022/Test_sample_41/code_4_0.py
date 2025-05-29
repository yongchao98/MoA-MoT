def find_rectangles(grid):
    lines = [line for line in grid.split('\n') if line.strip()]
    
    def is_horizontal_border(line):
        # Count consecutive '#' or '█' characters
        count = 0
        max_count = 0
        for c in line:
            if c in '#█':
                count += 1
            else:
                max_count = max(max_count, count)
                count = 0
        return max(max_count, count) >= 3  # Minimum width for a rectangle
    
    def find_matching_borders():
        rectangles = 0
        processed = set()
        
        # Find top borders
        for top in range(len(lines)-1):
            if is_horizontal_border(lines[top]):
                # Find left edge
                for left in range(len(lines[top])):
                    if lines[top][left] in '#█':
                        # Find right edge
                        for right in range(left + 2, len(lines[top])):
                            if lines[top][right] in '#█':
                                # Find bottom edge
                                for bottom in range(top + 1, len(lines)):
                                    if (is_horizontal_border(lines[bottom]) and
                                        all(lines[y][left] in '#█' and lines[y][right] in '#█' 
                                            for y in range(top, bottom + 1))):
                                        # Verify it's a new rectangle
                                        rect = (top, left, bottom, right)
                                        if rect not in processed:
                                            # Verify it has interior space
                                            has_space = False
                                            for y in range(top + 1, bottom):
                                                if ' ' in lines[y][left+1:right]:
                                                    has_space = True
                                                    break
                                            if has_space:
                                                processed.add(rect)
                                                rectangles += 1
                                        break
        return rectangles

    return find_matching_borders()

grid = """                                                                                
[... your grid here ...]
                                                                                """

print(find_rectangles(grid))
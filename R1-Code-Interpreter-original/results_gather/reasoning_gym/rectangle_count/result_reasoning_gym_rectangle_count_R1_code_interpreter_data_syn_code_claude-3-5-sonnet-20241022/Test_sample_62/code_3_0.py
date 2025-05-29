def find_rectangles(grid):
    # Convert grid to 2D array
    lines = grid.strip().splitlines()
    grid = [list(line) for line in lines]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0

    def find_top_left_corners():
        corners = []
        for i in range(height):
            for j in range(width):
                if grid[i][j] not in ['#', '█']:
                    continue
                    
                # Check if it's a top-left corner
                left_empty = j == 0 or grid[i][j-1] not in ['#', '█']
                top_empty = i == 0 or grid[i-1][j] not in ['#', '█']
                right_exists = j < width-1 and grid[i][j+1] in ['#', '█']
                bottom_exists = i < height-1 and grid[i+1][j] in ['#', '█']
                
                if left_empty and top_empty and right_exists and bottom_exists:
                    corners.append((i, j))
        return corners

    def trace_rectangle(start_i, start_j):
        # Find width
        width_count = 0
        j = start_j
        while j < width and grid[start_i][j] in ['#', '█']:
            width_count += 1
            j += 1

        # Find height
        height_count = 0
        i = start_i
        while i < height and grid[i][start_j] in ['#', '█']:
            height_count += 1
            i += 1

        # Verify rectangle
        end_i = start_i + height_count - 1
        end_j = start_j + width_count - 1
        
        # Check all corners exist
        if not (grid[start_i][end_j] in ['#', '█'] and
                grid[end_i][start_j] in ['#', '█'] and
                grid[end_i][end_j] in ['#', '█']):
            return None

        # Check all edges exist
        for i in range(start_i, end_i + 1):
            if grid[i][start_j] not in ['#', '█'] or grid[i][end_j] not in ['#', '█']:
                return None
        for j in range(start_j, end_j + 1):
            if grid[start_i][j] not in ['#', '█'] or grid[end_i][j] not in ['#', '█']:
                return None

        return (start_i, start_j, end_i, end_j)

    def is_new_rectangle(rect, found_rectangles):
        # Check if this rectangle is unique or part of overlapping rectangles
        for existing in found_rectangles:
            # If exactly the same, it's not new
            if rect == existing:
                return False
            
            # If it shares '█' characters but has different dimensions, it's new
            start_i, start_j, end_i, end_j = rect
            ex_start_i, ex_start_j, ex_end_i, ex_end_j = existing
            
            # Check for any '█' in the intersection
            has_overlap = False
            for i in range(max(start_i, ex_start_i), min(end_i + 1, ex_end_i + 1)):
                for j in range(max(start_j, ex_start_j), min(end_j + 1, ex_end_j + 1)):
                    if grid[i][j] == '█':
                        has_overlap = True
                        break
                if has_overlap:
                    break
            
            if has_overlap and (end_i - start_i != ex_end_i - ex_start_i or 
                              end_j - start_j != ex_end_j - ex_start_j):
                return True
        return True

    rectangles = []
    corners = find_top_left_corners()
    
    for corner in corners:
        rect = trace_rectangle(corner[0], corner[1])
        if rect and is_new_rectangle(rect, rectangles):
            rectangles.append(rect)

    return len(rectangles)

# Test grid
grid = """                                                                                
[... same grid as before ...]
                                       ##############################           """

print(find_rectangles(grid))
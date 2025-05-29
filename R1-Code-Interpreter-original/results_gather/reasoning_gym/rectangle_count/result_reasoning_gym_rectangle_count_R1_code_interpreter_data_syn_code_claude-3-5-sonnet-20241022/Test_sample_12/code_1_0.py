def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_rectangle(top, left, bottom, right):
        # Check corners
        if not all(grid[y][x] in '#█' for y, x in [
            (top, left), (top, right), 
            (bottom, left), (bottom, right)
        ]):
            return False
            
        # Check edges
        for x in range(left + 1, right):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        for y in range(top + 1, bottom):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
            
        # Check interior (should be empty except for possible overlaps)
        has_interior_hash = False
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if grid[y][x] == '#':
                    has_interior_hash = True
                    break
            if has_interior_hash:
                break
                
        return not has_interior_hash

    def find_next_corner(y, x):
        # Find the next corner point from current position
        while x < width and grid[y][x] not in '#█':
            x += 1
        return x if x < width else -1

    rectangles = []
    # Process each row
    for top in range(height):
        x = 0
        while x < width:
            left = find_next_corner(top, x)
            if left == -1:
                break
                
            # Find right edge
            right = left + 1
            while right < width and grid[top][right] in '#█':
                right += 1
            right -= 1
            
            # For each possible bottom
            for bottom in range(top + 1, height):
                if grid[bottom][left] in '#█' and grid[bottom][right] in '#█':
                    if is_valid_rectangle(top, left, bottom, right):
                        rectangles.append((top, left, bottom, right))
            
            x = right + 1

    # Remove duplicates
    unique_rectangles = []
    for rect in rectangles:
        is_unique = True
        for other in rectangles:
            if rect != other and (
                rect[0] >= other[0] and rect[1] >= other[1] and
                rect[2] <= other[2] and rect[3] <= other[3]
            ):
                # This rectangle is contained within another
                is_unique = False
                break
        if is_unique:
            unique_rectangles.append(rect)

    return len(unique_rectangles)

grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                         #########################              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       ############   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                    #####█#######################██#########█#  ",
    "                                    #    #                       ##         ##  ",
    "                                    #    #                       ##         ##  ",
    "                                    #    #                       ##         ##  ",
    "                                    #    #                       ##         ##  ",
    "                                    #    #                       ##         ##  ",
    "                                    #    #                       ##         ##  ",
    "                                    #    #                       ##         ##  ",
    "                                    #    #                       ##         ##  ",
    "                                    #####█#######################██#########█#  ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #                       ##         #   ",
    "                                         #      ##########       ##         #   ",
    "                                         #      #        #       ############   ",
    "                                         #      #        #       #              ",
    "                                         #      ##########       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #                       #              ",
    "                                         #########################              ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "            #######################                                             ",
    "            #                     #                                             ",
    "            #                     #                                             ",
    "            #                     #                                             ",
    "            #                     #                                             ",
    "            #                     #                                             ",
    "            #                     #                                             ",
    "            #               ######█###                                          ",
    "            #               #     #  #                                          ",
    "            #               ######█###                                          ",
    "            #                     #   ###########################               ",
    "            #                     #   #                         #               ",
    "            #                     #   #                         #               ",
    "            #######################   ###########################               ",
    "                                                                                "]

print(find_rectangles(grid))
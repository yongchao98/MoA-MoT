def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    height = len(grid)
    width = len(grid[0])
    rectangles = []

    # Function to check if a point is part of rectangle border
    def is_border(ch):
        return ch in '#█'

    # Function to validate if points form a valid rectangle
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True

    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if is_border(grid[top][left]):
                # Try all possible bottom-right combinations
                for bottom in range(top, height):
                    for right in range(left, width):
                        if (is_border(grid[top][right]) and 
                            is_border(grid[bottom][left]) and 
                            is_border(grid[bottom][right])):
                            if is_valid_rectangle(top, left, bottom, right):
                                rectangles.append((top, left, bottom, right))

    # Count overlapping rectangles
    overlaps = set()
    for i, rect1 in enumerate(rectangles):
        for j, rect2 in enumerate(rectangles[i+1:], i+1):
            top1, left1, bottom1, right1 = rect1
            top2, left2, bottom2, right2 = rect2
            
            # Check for overlap
            if (left1 <= right2 and right1 >= left2 and 
                top1 <= bottom2 and bottom1 >= top2):
                # Check if there's a '█' character in the overlap
                for y in range(max(top1, top2), min(bottom1, bottom2) + 1):
                    for x in range(max(left1, left2), min(right1, right2) + 1):
                        if grid[y][x] == '█':
                            overlaps.add((i, j))
                            break
                    if (i, j) in overlaps:
                        break

    # Total rectangles = unique rectangles + overlapping pairs
    total_rectangles = len(rectangles) - len(overlaps)
    print(total_rectangles)

# Create the grid string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                          ###########################           
                                          #                         #           
                                          #                   ######█########   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #                   #     #       #   
                                          #             ######█#####█########   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     #     #      ##   
                                          #             #     ######█######█#   
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             #           #      #    
                                          #             ############█#######    
                                          #                         #           
                                          #                         #           
                                          #                         #           
                                          ###########################           
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)
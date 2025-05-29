def find_rectangles(grid):
    # Convert string grid to 2D array
    grid = [list(row) for row in grid.splitlines()]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = set()

    def verify_rectangle(top, left, bottom, right):
        # Verify the corners first
        if grid[top][left] not in ['#', '█'] or grid[top][right] not in ['#', '█'] or \
           grid[bottom][left] not in ['#', '█'] or grid[bottom][right] not in ['#', '█']:
            return False

        # Check that interior is empty (except for borders)
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if grid[y][x] not in [' ', '#', '█']:
                    return False
                
        # Verify top and bottom edges are complete
        for x in range(left + 1, right):
            if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return False

        # Verify left and right edges are complete
        for y in range(top + 1, bottom):
            if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return False

        # Check that it's not just a line
        if bottom - top < 1 or right - left < 1:
            return False

        return True

    def find_complete_rectangle(y, x):
        if grid[y][x] not in ['#', '█']:
            return None
        
        # Find potential right edge
        right = x + 1
        while right < width and grid[y][right] in ['#', '█']:
            right += 1
        right -= 1

        # Find potential bottom edge
        bottom = y + 1
        while bottom < height and grid[bottom][x] in ['#', '█']:
            bottom += 1
        bottom -= 1

        if verify_rectangle(y, x, bottom, right):
            return (y, x, bottom, right)
        return None

    # Scan for rectangles
    for y in range(height):
        for x in range(width):
            rect = find_complete_rectangle(y, x)
            if rect:
                # Only add if it's a top-left corner
                if (y == 0 or grid[y-1][x] not in ['#', '█']) and \
                   (x == 0 or grid[y][x-1] not in ['#', '█']):
                    rectangles.add(rect)

    return len(rectangles)

# Test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                               #####            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                         ###############       #   #            
                                         #             #       #   #            
                                         #             #       #   #            
                                         #             #       #   #            
                                         #             #       #   #            
                                         #             #       #####            
                                         #             #       ## ##            
                                         #             #       ## ##            
                                         #             #       #####            
                                         #             #       #   #            
                                         #             #       #   #            
                                         #             #       #   #            
                                         #             #       #   #            
                                         #             #       #   #            
                                         ###############       #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #   #            
                                                               #####            
                                                                                
              ########################################                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #                                      #                          
              #               ############           #                          
              #               #          #           #                          
              #               #          #           #                          
              #               #          #           #                          
              #               #          #           #                          
              #               #          #           #                          
              #               #          #           #                          
              #               ############           #                          
              ########################################                          """

print(find_rectangles(grid))
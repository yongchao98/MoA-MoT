def find_rectangles(grid):
    # Convert string grid to 2D array
    grid = [list(row) for row in grid.splitlines()]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = []

    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return False
        return True

    # Scan for rectangles
    for y in range(height):
        for x in range(width):
            # Find potential top-left corners
            if grid[y][x] in ['#', '█']:
                # Look for bottom-right corner
                right = x + 1
                while right < width and grid[y][right] in ['#', '█']:
                    right += 1
                right -= 1

                bottom = y + 1
                while bottom < height and grid[bottom][x] in ['#', '█']:
                    bottom += 1
                bottom -= 1

                # Verify if it forms a valid rectangle
                if is_valid_rectangle(y, x, bottom, right):
                    rectangles.append((y, x, bottom, right))

    # Count unique rectangles
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
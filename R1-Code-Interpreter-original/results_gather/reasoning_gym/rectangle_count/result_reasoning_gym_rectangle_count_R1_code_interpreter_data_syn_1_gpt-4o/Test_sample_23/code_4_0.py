# Define the grid as a list of strings
grid = [
    "                       ######################                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       # ############       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # ############       #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       ######################                                   "
]

def count_rectangles(grid):
    # Initialize a set to store unique rectangles
    rectangles = set()
    
    # Get the dimensions of the grid
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    # Create a visited matrix
    visited = [[False] * width for _ in range(height)]
    
    # Iterate over each cell in the grid
    for y in range(height):
        for x in range(width):
            # Check for the top-left corner of a rectangle
            if grid[y][x] == '#' and not visited[y][x]:
                # Find the bottom-right corner of the rectangle
                x2, y2 = x, y
                while x2 < width and grid[y][x2] == '#':
                    x2 += 1
                while y2 < height and grid[y2][x] == '#':
                    y2 += 1
                
                # Check if the rectangle is valid by ensuring the bottom-right corner is also '#'
                if grid[y2-1][x2-1] == '#':
                    # Mark the rectangle as visited
                    for i in range(y, y2):
                        for j in range(x, x2):
                            visited[i][j] = True
                    
                    # Add the rectangle to the set
                    rectangles.add((x, y, x2 - 1, y2 - 1))
    
    # Return the number of unique rectangles
    return len(rectangles)

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)
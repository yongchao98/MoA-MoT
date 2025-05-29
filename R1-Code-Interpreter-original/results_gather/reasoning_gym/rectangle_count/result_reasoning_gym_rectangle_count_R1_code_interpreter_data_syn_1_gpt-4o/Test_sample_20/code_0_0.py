def count_rectangles(grid):
    # Function to count rectangles in the grid
    def is_rectangle(r, c, visited):
        # Check if the current position is the start of a rectangle
        if grid[r][c] != '#' and grid[r][c] != '█':
            return False
        # Check if this position is already visited
        if visited[r][c]:
            return False
        # Mark the current position as visited
        visited[r][c] = True
        # Determine the size of the rectangle
        width, height = 0, 0
        # Calculate width
        while c + width < len(grid[0]) and (grid[r][c + width] == '#' or grid[r][c + width] == '█'):
            width += 1
        # Calculate height
        while r + height < len(grid) and (grid[r + height][c] == '#' or grid[r + height][c] == '█'):
            height += 1
        # Mark the entire rectangle as visited
        for i in range(height):
            for j in range(width):
                visited[r + i][c + j] = True
        return True

    visited = [[False] * len(grid[0]) for _ in range(len(grid))]
    rectangle_count = 0

    for r in range(len(grid)):
        for c in range(len(grid[0])):
            if is_rectangle(r, c, visited):
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    " " * 48 + "####" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "#  #" + " " * 24,
    " " * 48 + "####" + " " * 24,
    " " * 52 + "###" + " " * 20,
    " " * 52 + "# #" + " " * 20,
    " " * 52 + "# #" + " " * 20,
    " " * 52 + "# #" + " " * 20,
    " " * 52 + "# #" + " " * 20,
    " " * 52 + "# #" + " " * 20,
    " " * 52 + "###" + " " * 20,
    " " * 48 + "########" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 48 + "#      #" + " " * 16,
    " " * 6 + "##############################################█######█##############" + " " * 6,
    " " * 6 + "#                                             #      #             #" + " " * 6,
    " " * 6 + "#                                             #      #             #" + " " * 6,
    " " * 6 + "#                                             #      #             #" + " " * 6,
    " " * 6 + "#                                             #      #             #" + " " * 6,
    " " * 6 + "#                                             #      #             #" + " " * 6,
    " " * 6 + "#                                             #      #             #" + " " * 6,
    " " * 6 + "#                                             #      #             #" + " " * 6,
    " " * 6 + "######█                                             #      #             #" + " " * 6,
    " " * 6 + "#     █                                             #      #             #" + " " * 6,
    " " * 6 + "#     █                                             #      #             #" + " " * 6,
    " " * 6 + "#     █                                             #      #             #" + " " * 6,
    " " * 6 + "#     █                                             #      #             #" + " " * 6,
    " " * 6 + "#     █#############################################████████##############" + " " * 6,
    " " * 6 + "#     #" + " " * 64,
    " " * 6 + "#     #" + " " * 64,
    " " * 6 + "#######" + " " * 64,
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)
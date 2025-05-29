def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def is_valid(x, y):
        return 0 <= x < rows and 0 <= y < cols
    
    def flood_fill_check(x, y, visited, background):
        if not is_valid(x, y) or (x, y) in visited:
            return True
        if grid[x][y] == 8:
            return True
        if x == 0 or x == rows-1 or y == 0 or y == cols-1:
            return False
        visited.add((x, y))
        
        # Check all 4 directions
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            if not flood_fill_check(x + dx, y + dy, visited, background):
                return False
        return True
    
    def flood_fill_apply(x, y, background):
        if not is_valid(x, y) or grid[x][y] == 8 or output[x][y] == 0:
            return
        if grid[x][y] != background:
            return
            
        output[x][y] = 0
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            flood_fill_apply(x + dx, y + dy, background)
    
    # Check each cell
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 8 and output[i][j] != 0:
                visited = set()
                background = grid[i][j]
                if flood_fill_check(i, j, visited, background):
                    flood_fill_apply(i, j, background)
    
    # Convert to string format
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = [
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,8,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,5,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,5,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,8,8,8,5,5,5,5],
    [5,8,8,8,8,5,5,5,5,5,5,5,5,5,5],
    [5,8,5,5,8,5,5,5,8,8,8,8,5,5,5],
    [5,8,5,5,8,5,5,5,8,5,5,8,5,5,5],
    [5,8,5,8,8,5,5,5,8,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,8,8,8,8,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
]

find_enclosed_areas(test_input)
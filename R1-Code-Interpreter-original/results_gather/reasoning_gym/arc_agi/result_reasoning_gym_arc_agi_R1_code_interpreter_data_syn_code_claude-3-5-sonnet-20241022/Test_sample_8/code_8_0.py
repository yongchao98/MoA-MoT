def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def flood_fill(i, j, visited):
        if (i < 0 or i >= rows or j < 0 or j >= cols or 
            (i,j) in visited or grid[i][j] == 8):
            return
        
        visited.add((i,j))
        
        # Check if this position should be filled
        left_8 = False
        right_8 = False
        
        # Look left for 8
        for y in range(j-1, -1, -1):
            if grid[i][y] == 8:
                left_8 = True
                break
            if grid[i][y] != grid[i][j]:
                break
                
        # Look right for 8
        for y in range(j+1, cols):
            if grid[i][y] == 8:
                right_8 = True
                break
            if grid[i][y] != grid[i][j]:
                break
        
        # If bounded horizontally, check for vertical connection
        if left_8 and right_8:
            has_vertical = False
            # Check up
            for x in range(i-1, -1, -1):
                if grid[x][j] == 8:
                    has_vertical = True
                    break
                if grid[x][j] != grid[i][j]:
                    break
            # Check down
            if not has_vertical:
                for x in range(i+1, rows):
                    if grid[x][j] == 8:
                        has_vertical = True
                        break
                    if grid[x][j] != grid[i][j]:
                        break
            
            if has_vertical:
                output[i][j] = 0
                # Continue flood fill in all directions
                for ni, nj in [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]:
                    flood_fill(ni, nj, visited)
    
    # Process each cell
    visited = set()
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 8 and (i,j) not in visited:
                flood_fill(i, j, visited)
    
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
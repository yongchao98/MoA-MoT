def find_continuous_block(grid):
    rows = len(grid)
    cols = len(grid[0])
    visited = [[False] * cols for _ in range(rows)]
    
    def is_valid(x, y):
        return 0 <= x < rows and 0 <= y < cols
    
    def dfs(x, y, val, block):
        if not is_valid(x, y) or visited[x][y] or grid[x][y] not in [4, 5]:
            return
        
        visited[x][y] = True
        block.append((x, y))
        
        # Check all 4 directions
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            new_x, new_y = x + dy, y + dx
            dfs(new_x, new_y, val, block)
    
    # Find the first significant block
    best_block = []
    for i in range(rows):
        for j in range(cols):
            if not visited[i][j] and grid[i][j] in [4, 5]:
                current_block = []
                dfs(i, j, grid[i][j], current_block)
                if len(current_block) >= 4 and (not best_block or len(current_block) > len(best_block)):
                    best_block = current_block
    
    return best_block

def solve(input_grid):
    # Convert input string to 2D list
    grid = [[int(num) for num in row.split()] for row in input_grid.strip().split('\n')]
    
    # Create output grid as a copy of input grid
    output_grid = [row[:] for row in grid]
    
    # Find the first significant continuous block
    block = find_continuous_block(grid)
    
    # Replace the block with 7s if found
    if block:
        min_x = min(x for x, y in block)
        max_x = max(x for x, y in block)
        min_y = min(y for x, y in block)
        max_y = max(y for x, y in block)
        
        # Check if block is significant (at least 2x2)
        if max_x - min_x >= 1 and max_y - min_y >= 1:
            for x in range(min_x, max_x + 1):
                for y in range(min_y, max_y + 1):
                    output_grid[x][y] = 7
    
    # Convert back to string format
    return '\n'.join(' '.join(str(num) for num in row) for row in output_grid)

# Test input
test_input = """4 4 5 5
5 4 4 4
5 4 4 5
5 4 4 4
4 5 4 4
4 5 4 5
5 4 5 5
4 4 4 5
4 5 4 5
5 5 5 4
4 5 4 5
4 4 4 5
4 5 5 5
4 4 5 5
4 4 5 5
5 5 5 4
5 4 4 5
5 4 4 4
5 4 4 5
5 5 4 4
4 5 5 5
4 5 5 5
5 5 5 5
5 5 4 5"""

print(solve(test_input))
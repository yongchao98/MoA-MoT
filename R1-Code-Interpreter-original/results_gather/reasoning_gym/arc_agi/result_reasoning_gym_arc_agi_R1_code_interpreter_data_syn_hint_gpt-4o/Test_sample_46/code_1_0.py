def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[7] * cols for _ in range(rows)]
    visited = [[False] * cols for _ in range(rows)]
    
    def flood_fill(r, c, num):
        # Use a stack for iterative flood fill
        stack = [(r, c)]
        block = []
        while stack:
            x, y = stack.pop()
            if 0 <= x < rows and 0 <= y < cols and not visited[x][y] and input_grid[x][y] == num:
                visited[x][y] = True
                block.append((x, y))
                # Add neighboring cells to the stack
                stack.extend([(x-1, y), (x+1, y), (x, y-1), (x, y+1)])
        return block
    
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != 7 and not visited[r][c]:
                block = flood_fill(r, c, input_grid[r][c])
                if len(block) > 1:  # Retain only if it's a contiguous block
                    for x, y in block:
                        output_grid[x][y] = input_grid[x][y]
    
    return output_grid

# Test input grid
input_grid = [
    [6, 7, 7, 7, 6, 7, 7, 7, 7, 6],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 6, 7, 6, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [6, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 6, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 6, 7],
    [7, 7, 6, 7, 7, 6, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
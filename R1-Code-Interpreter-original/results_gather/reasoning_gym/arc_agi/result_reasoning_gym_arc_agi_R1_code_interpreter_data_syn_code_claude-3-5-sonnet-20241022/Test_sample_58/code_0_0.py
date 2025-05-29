def find_bounds(grid):
    rows, cols = len(grid), len(grid[0])
    min_row, max_row, min_col, max_col = rows, 0, cols, 0
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 0:
                min_row = min(min_row, i)
                max_row = max(max_row, i)
                min_col = min(min_col, j)
                max_col = max(max_col, j)
    return min_row, max_row, min_col, max_col

def find_components(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    components = []
    
    def dfs(r, c, val):
        if (r, c) in visited or r < 0 or r >= rows or c < 0 or c >= cols or grid[r][c] != val:
            return []
        visited.add((r, c))
        points = [(r, c)]
        for dr, dc in [(0,1), (1,0), (0,-1), (-1,0), (1,1), (-1,-1), (1,-1), (-1,1)]:
            points.extend(dfs(r+dr, c+dc, val))
        return points
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 0 and (i,j) not in visited:
                comp = dfs(i, j, grid[i][j])
                if comp:
                    components.append((grid[i][j], comp))
    return components

def create_output_grid(input_grid):
    # Find bounds of non-zero elements
    min_r, max_r, min_c, max_c = find_bounds(input_grid)
    components = find_components(input_grid)
    
    # Find corners with 8s
    corners = []
    rows, cols = len(input_grid), len(input_grid[0])
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 8:
                corners.append((i,j))
    
    # Calculate output dimensions
    out_height = max(abs(max_r - min_r) // 2 + 3, len(corners) // 2 + 1)
    out_width = max(abs(max_c - min_c) // 2 + 3, len(corners) // 2 + 1)
    
    # Create output grid
    output = [[0 for _ in range(out_width)] for _ in range(out_height)]
    
    # Place corners
    if corners:
        output[0][0] = 8
        output[0][-1] = 8
        output[-1][0] = 8
        output[-1][-1] = 8
    
    # Place components in compressed form
    for val, points in components:
        if val != 8:
            min_x = min(x for x,y in points)
            max_x = max(x for x,y in points)
            min_y = min(y for x,y in points)
            max_y = max(y for x,y in points)
            
            # Calculate center position in output grid
            center_r = (out_height - 2) // 2
            center_c = (out_width - 2) // 2
            
            block_height = (max_x - min_x + 2) // 2
            block_width = (max_y - min_y + 2) // 2
            
            # Place the block
            for i in range(block_height):
                for j in range(block_width):
                    r = center_r - block_height//2 + i
                    c = center_c - block_width//2 + j
                    if 0 <= r < out_height and 0 <= c < out_width:
                        output[r][c] = val
    
    # Format output
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,7,7,7,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,4,0,2,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8]
]

create_output_grid(test_input)
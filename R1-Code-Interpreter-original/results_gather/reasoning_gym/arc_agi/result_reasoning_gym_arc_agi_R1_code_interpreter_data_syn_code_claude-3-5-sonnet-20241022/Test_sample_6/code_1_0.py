def is_valid_3x3(grid, r, c):
    rows, cols = len(grid), len(grid[0])
    if r < 1 or r >= rows-1 or c < 1 or c >= cols-1:
        return False
    if grid[r][c] != 9:
        return False
    # Check all surrounding positions must be 6
    return all(grid[r+dr][c+dc] == 6 
              for dr in [-1,0,1] 
              for dc in [-1,0,1] 
              if (dr,dc) != (0,0))

def find_component(grid, r, c, visited):
    rows, cols = len(grid), len(grid[0])
    if (r,c) in visited or grid[r][c] not in [6,9]:
        return set()
    
    component = {(r,c)}
    visited.add((r,c))
    stack = [(r,c)]
    
    while stack:
        curr_r, curr_c = stack.pop()
        for dr in [-1,0,1]:
            for dc in [-1,0,1]:
                new_r, new_c = curr_r + dr, curr_c + dc
                if (0 <= new_r < rows and 
                    0 <= new_c < cols and 
                    (new_r,new_c) not in visited and 
                    grid[new_r][new_c] in [6,9]):
                    component.add((new_r,new_c))
                    visited.add((new_r,new_c))
                    stack.append((new_r,new_c))
    return component

def transform_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = [[8]*cols for _ in range(rows)]
    visited = set()
    
    # Find all components starting from bottom row
    for j in range(cols):
        if grid[rows-1][j] in [6,9] and (rows-1,j) not in visited:
            # Get the component
            component = find_component(grid, rows-1, j, visited)
            
            # Check if component has a valid 3x3 pattern
            has_valid_pattern = any(is_valid_3x3(grid, r, c) 
                                  for r,c in component)
            
            # If component has valid pattern, keep it in output
            if has_valid_pattern:
                for r,c in component:
                    output[r][c] = grid[r][c]
    
    return output

# Test input
test_input = [
    [8,8,8,8,8,8,8,8,6,6,6,8],
    [8,6,6,6,6,6,6,8,6,6,6,8],
    [8,6,9,6,6,6,6,8,8,8,8,8],
    [8,6,6,6,6,9,6,8,6,6,6,6],
    [8,6,6,6,6,6,6,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,6,9,6],
    [8,8,8,8,6,6,6,8,6,6,6,6],
    [8,8,8,8,6,6,9,8,8,8,8,8],
    [8,8,8,8,6,6,6,8,8,8,8,8],
    [6,6,6,8,9,6,6,8,8,8,8,8],
    [9,6,6,8,8,8,8,8,6,9,6,9],
    [6,6,6,8,8,8,8,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,9,6,6]
]

result = transform_grid(test_input)
print("Final output:")
for row in result:
    print(' '.join(map(str, row)))
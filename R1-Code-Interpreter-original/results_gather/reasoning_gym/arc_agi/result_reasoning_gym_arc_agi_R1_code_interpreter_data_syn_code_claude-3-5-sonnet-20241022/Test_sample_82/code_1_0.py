def find_connected_components(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    components = []
    
    def dfs(r, c, component):
        if (r < 0 or r >= rows or c < 0 or c >= cols or 
            (r,c) in visited or grid[r][c] == 1):
            return
        visited.add((r,c))
        component.add((r,c))
        val = grid[r][c]
        # Check all 8 directions (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                new_r, new_c = r + dr, c + dc
                if (new_r >= 0 and new_r < rows and new_c >= 0 and new_c < cols and 
                    grid[new_r][new_c] != 1 and (new_r, new_c) not in visited):
                    dfs(new_r, new_c, component)
    
    for i in range(rows):
        for j in range(cols):
            if (i,j) not in visited and grid[i][j] != 1:
                current_component = set()
                dfs(i, j, current_component)
                if current_component:
                    components.append(current_component)
    
    return components

def get_output_format(grid):
    components = find_connected_components(grid)
    output_lines = []
    
    for component in components:
        numbers = sorted(set(grid[r][c] for r,c in component))
        output_lines.append(numbers)
    
    # Find the size of the largest component
    max_size = max(len(component) for component in components)
    
    # Format output
    result = []
    for numbers in output_lines:
        result.append(' '.join(map(str, numbers)))
    
    # Return unique lines and repetition count
    unique_lines = sorted(set(result))
    return unique_lines, max_size

# Process test input
test_input = """1 8 1 1 1 1 1 1 1 1
8 8 8 1 1 1 1 6 1 1
1 1 1 1 1 1 1 6 6 1
1 1 1 1 1 1 1 6 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 3 3 1 1 1 1 1
1 1 1 1 1 1 1 1 1 0
1 1 1 1 1 9 9 1 1 0
1 1 1 1 1 1 9 1 1 0
1 1 1 1 1 1 9 1 1 1"""

grid = [[int(x) for x in row.split()] for row in test_input.strip().split('\n')]
lines, repeat_count = get_output_format(grid)
for line in lines:
    print(line)
print(f"Repeat count: {repeat_count}")
def find_left_connected_components(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    components = []
    
    def dfs(r, c, component):
        if (r < 0 or r >= rows or c < 0 or c >= cols or 
            (r,c) in visited or grid[r][c] == 1):
            return
        visited.add((r,c))
        component.add((r,c))
        # Check all 8 directions (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                new_r, new_c = r + dr, c + dc
                if (new_r >= 0 and new_r < rows and new_c >= 0 and new_c < cols and 
                    grid[new_r][new_c] != 1 and (new_r, new_c) not in visited):
                    dfs(new_r, new_c, component)

    # First check leftmost column
    left_components = []
    for i in range(rows):
        if grid[i][0] != 1 and (i,0) not in visited:
            current_component = set()
            dfs(i, 0, current_component)
            if current_component:
                numbers = sorted(set(grid[r][c] for r,c in current_component))
                left_components.append(numbers)
    
    return left_components

def process_grid(input_grid):
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    components = find_left_connected_components(grid)
    
    # Format output lines
    output_lines = []
    for numbers in components:
        output_lines.append(' '.join(map(str, numbers)))
    
    return output_lines

# Test input
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

result = process_grid(test_input)
for line in result:
    print(line)
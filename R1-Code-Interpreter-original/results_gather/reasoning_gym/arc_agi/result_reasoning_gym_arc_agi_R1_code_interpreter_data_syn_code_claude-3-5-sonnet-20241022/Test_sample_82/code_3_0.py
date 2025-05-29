def find_connected_groups(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    groups = []
    
    def dfs(r, c, group):
        if (r < 0 or r >= rows or c < 0 or c >= cols or 
            (r,c) in visited or grid[r][c] == 1):
            return
        visited.add((r,c))
        group.add((r,c))
        val = grid[r][c]
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                new_r, new_c = r + dr, c + dc
                if (new_r >= 0 and new_r < rows and new_c >= 0 and new_c < cols and 
                    grid[new_r][new_c] != 1 and (new_r, new_c) not in visited):
                    dfs(new_r, new_c, group)

    # First check leftmost column
    left_edge_groups = []
    for i in range(rows):
        if grid[i][0] != 1 and (i,0) not in visited:
            current_group = set()
            dfs(i, 0, current_group)
            if current_group:
                left_edge_groups.append(current_group)
    
    # If no groups from left edge, find all connected groups
    if not left_edge_groups:
        for i in range(rows):
            for j in range(cols):
                if grid[i][j] != 1 and (i,j) not in visited:
                    current_group = set()
                    dfs(i, j, current_group)
                    if current_group:
                        groups.append(current_group)
        return groups
    return left_edge_groups

def process_grid(input_grid):
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    groups = find_connected_groups(grid)
    
    # Get numbers from groups
    group_numbers = []
    max_size = 0
    for group in groups:
        numbers = sorted(set(grid[r][c] for r,c in group))
        group_numbers.append(numbers)
        max_size = max(max_size, len(group))
    
    # Format output line
    output_line = ' '.join(map(str, sorted(set([num for group in group_numbers for num in group]))))
    
    return [output_line] * max_size

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
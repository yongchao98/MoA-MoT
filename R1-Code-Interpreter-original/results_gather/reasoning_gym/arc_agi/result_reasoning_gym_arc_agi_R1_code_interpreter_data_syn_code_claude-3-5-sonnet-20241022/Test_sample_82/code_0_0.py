def find_connected_groups(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    groups = []
    
    def dfs(r, c, num, group):
        if (r < 0 or r >= rows or c < 0 or c >= cols or 
            (r,c) in visited or grid[r][c] != num):
            return
        visited.add((r,c))
        group.add((r,c))
        for dr, dc in [(0,1), (1,0), (0,-1), (-1,0), (1,1), (-1,-1), (1,-1), (-1,1)]:
            dfs(r+dr, c+dc, num, group)
    
    for i in range(rows):
        for j in range(cols):
            if (i,j) not in visited and grid[i][j] != 1:
                current_group = set()
                dfs(i, j, grid[i][j], current_group)
                if current_group:
                    groups.append(current_group)
    
    return groups

def process_grid(input_grid):
    # Convert input string to 2D grid
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Find connected groups
    groups = find_connected_groups(grid)
    
    # Extract unique numbers from each group
    numbers = []
    for group in groups:
        nums = sorted(set(grid[r][c] for r,c in group))
        if nums:
            numbers.extend(nums)
    
    # Remove duplicates and sort
    numbers = sorted(set(numbers))
    
    # Determine number of repetitions (using length of first example as default)
    repetitions = max(len(groups), 3)
    
    # Generate output
    output = [' '.join(map(str, numbers)) for _ in range(repetitions)]
    return output

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
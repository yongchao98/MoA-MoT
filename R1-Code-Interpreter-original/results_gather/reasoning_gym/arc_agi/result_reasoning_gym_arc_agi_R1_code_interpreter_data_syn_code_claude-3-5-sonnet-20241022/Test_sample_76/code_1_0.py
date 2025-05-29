def analyze_grid(input_grid):
    # Convert input string to 2D grid
    rows = input_grid.strip().split('\n')
    grid = [[int(x) for x in row.split()] for row in rows]
    
    # Find special numbers (non-3)
    special_nums = set()
    for row in grid:
        for num in row:
            if num != 3:
                special_nums.add(num)
    
    # For each special number, find its largest cluster
    clusters = {}
    rows, cols = len(grid), len(grid[0])
    
    def count_connected(r, c, num, visited):
        if (r < 0 or r >= rows or c < 0 or c >= cols or 
            grid[r][c] != num or (r,c) in visited):
            return 0
        visited.add((r,c))
        count = 1
        # Check all 8 directions
        for dr in [-1,0,1]:
            for dc in [-1,0,1]:
                if dr == 0 and dc == 0:
                    continue
                count += count_connected(r+dr, c+dc, num, visited)
        return count
    
    # Find largest cluster for each special number
    for num in special_nums:
        max_cluster = 0
        for i in range(rows):
            for j in range(cols):
                if grid[i][j] == num:
                    cluster_size = count_connected(i, j, num, set())
                    max_cluster = max(max_cluster, cluster_size)
        clusters[num] = max_cluster
    
    # Analyze distribution in 3x3 regions
    height, width = rows//3, cols//3
    output = []
    
    for i in range(3):
        row = []
        for j in range(3):
            region_nums = set()
            for r in range(i*height, (i+1)*height):
                for c in range(j*width, (j+1)*width):
                    if grid[r][c] != 3:
                        region_nums.add(grid[r][c])
            
            # If region has special numbers, use most significant one
            if region_nums:
                # Prefer numbers that form larger clusters
                significant_num = max(region_nums, key=lambda x: clusters[x])
                row.append(significant_num)
            else:
                row.append(3)
        output.append(row)
    
    # Print the output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = """3 3 3 3 3 2 3 3 3 3 2 2 3 3 2 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 2 3 3 3 3 3 3 3 3
2 3 3 3 2 3 3 3 3 2 3 3 2 3 3 3 3 3 3 3 3 2 3 3 3 3
3 3 3 3 3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 2 3 3 3 3 2 3 2 3 3 3 3 2 2 3 2 3
2 3 3 3 3 3 2 3 2 3 3 3 4 4 4 4 3 3 3 3 2 3 3 2 3 3
3 3 3 3 3 3 2 3 3 3 3 3 4 4 4 4 3 2 3 3 3 3 3 3 3 2
2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 2 3 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 2 3 3 3 3 2 3 3 2 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 3 3 2 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3
3 2 3 3 3 3 3 3 3 2 2 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3
3 3 2 3 2 3 3 3 3 3 3 2 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 2 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 3 3 3 2 3 3 3 3 3 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 2 3 3 3 3 2 3 3 3 3
3 3 3 3 2 2 2 2 3 3 3 3 2 3 3 3 3 2 3 3 3 3 3 3 3 3
3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 3 3 2 3 3 3
3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3"""

analyze_grid(test_input)
def find_pairs_by_region(grid):
    rows, cols = len(grid), len(grid[0])
    pairs = []
    visited = set()
    
    # Helper function to check if position is valid
    def is_valid(r, c):
        return 0 <= r < rows and 0 <= c < cols
    
    # Helper function to find pair
    def find_pair(r, c):
        if (r, c) in visited or grid[r][c] == 5:
            return None
        num = grid[r][c]
        # Check right and down
        for dr, dc in [(0,1), (1,0)]:
            nr, nc = r + dr, c + dc
            if is_valid(nr, nc) and grid[nr][nc] == num:
                visited.add((r, c))
                visited.add((nr, nc))
                return {
                    'num': num,
                    'type': 'horizontal' if dr == 0 else 'vertical',
                    'row': r,
                    'col': c,
                    'region': 'top' if r < rows//3 else ('middle' if r < 2*rows//3 else 'bottom')
                }
        return None

    # Scan grid for pairs
    for r in range(rows):
        for c in range(cols):
            pair = find_pair(r, c)
            if pair:
                pairs.append(pair)
    
    return pairs

def solve_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Find all pairs and categorize by region
    pairs = find_pairs_by_region(grid)
    
    # Sort pairs by row position within each region
    pairs.sort(key=lambda x: (x['row'], x['col']))
    
    # Create 4x4 output grid
    output = [[5]*4 for _ in range(4)]
    
    # Process pairs and place them in output grid
    row_idx = 0
    for pair in pairs:
        if pair['col'] > len(grid[0])//2:  # Right side pairs
            if row_idx < 2:  # Top two rows
                output[row_idx][2] = pair['num']
                output[row_idx][3] = 5 if pair['type'] == 'horizontal' else pair['num']
        else:  # Left side pairs
            if row_idx < 4:  # All rows
                output[row_idx][0] = pair['num']
                output[row_idx][1] = 5 if pair['type'] == 'horizontal' else pair['num']
        row_idx = (row_idx + 1) if pair['col'] <= len(grid[0])//2 else row_idx
    
    # Format and print output
    result = '\n'.join(' '.join(str(x) for x in row) for row in output)
    print(result)

# Test input
test_input = """5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 7 7 5 5
5 5 5 5 5 5 5 5 5 5 7 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 2 5 5 5 5 5 5 5
5 5 5 5 5 2 2 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 0 5 5 5 5 5 1 1 5 5
5 5 0 0 5 5 5 5 5 1 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5"""

solve_grid(test_input)
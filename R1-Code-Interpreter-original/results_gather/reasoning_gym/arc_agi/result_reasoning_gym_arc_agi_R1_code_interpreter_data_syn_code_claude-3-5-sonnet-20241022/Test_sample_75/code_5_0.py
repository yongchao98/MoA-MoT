def find_pairs_by_region(grid):
    rows, cols = len(grid), len(grid[0])
    pairs = []
    visited = set()
    
    def find_pair(r, c):
        if (r, c) in visited or grid[r][c] == 5:
            return None
        num = grid[r][c]
        # Check right and down
        for dr, dc in [(0,1), (1,0)]:
            nr, nc = r + dr, c + dc
            if (0 <= nr < rows and 0 <= nc < cols and 
                grid[nr][nc] == num and (nr, nc) not in visited):
                visited.add((r, c))
                visited.add((nr, nc))
                return {
                    'num': num,
                    'type': 'horizontal' if dr == 0 else 'vertical',
                    'row': r,
                    'col': c
                }
        return None

    # Find all pairs
    for r in range(rows):
        for c in range(cols):
            pair = find_pair(r, c)
            if pair:
                pairs.append(pair)
    
    # Separate pairs into left and right
    left_pairs = [p for p in pairs if p['col'] < cols//2]
    right_pairs = [p for p in pairs if p['col'] >= cols//2]
    
    # Sort each group by row position
    left_pairs.sort(key=lambda x: x['row'])
    right_pairs.sort(key=lambda x: x['row'])
    
    return left_pairs, right_pairs

def solve_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Get sorted pairs
    left_pairs, right_pairs = find_pairs_by_region(grid)
    
    # Create 4x4 output grid
    output = [[5]*4 for _ in range(4)]
    
    # Place left pairs
    if len(left_pairs) >= 1:
        # First left pair
        output[0][0] = left_pairs[0]['num']
        output[0][1] = 5 if left_pairs[0]['type'] == 'horizontal' else left_pairs[0]['num']
    if len(left_pairs) >= 2:
        # Second left pair
        output[2][0] = left_pairs[1]['num']
        output[2][1] = 5 if left_pairs[1]['type'] == 'horizontal' else left_pairs[1]['num']
        if left_pairs[1]['type'] == 'vertical':
            output[3][0] = left_pairs[1]['num']
            output[3][1] = 5
    
    # Place right pairs
    if len(right_pairs) >= 1:
        # First right pair
        output[0][2] = right_pairs[0]['num']
        output[0][3] = 5 if right_pairs[0]['type'] == 'horizontal' else right_pairs[0]['num']
        if right_pairs[0]['type'] == 'vertical':
            output[1][2] = right_pairs[0]['num']
            output[1][3] = 5
    if len(right_pairs) >= 2:
        # Second right pair
        output[2][2] = right_pairs[1]['num']
        output[2][3] = 5 if right_pairs[1]['type'] == 'horizontal' else right_pairs[1]['num']
        if right_pairs[1]['type'] == 'vertical':
            output[3][2] = right_pairs[1]['num']
            output[3][3] = 5
    
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
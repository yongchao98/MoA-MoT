def find_and_sort_pairs(grid):
    rows, cols = len(grid), len(grid[0])
    pairs = []
    visited = set()
    
    # Find all pairs
    for r in range(rows):
        for c in range(cols):
            if (r, c) not in visited and grid[r][c] != 5:
                num = grid[r][c]
                # Check right and down
                for dr, dc in [(0,1), (1,0)]:
                    nr, nc = r + dr, c + dc
                    if (0 <= nr < rows and 0 <= nc < cols and 
                        grid[nr][nc] == num and (nr, nc) not in visited):
                        pairs.append({
                            'num': num,
                            'row': r,
                            'col': c,
                            'type': 'horizontal' if dr == 0 else 'vertical'
                        })
                        visited.add((r, c))
                        visited.add((nr, nc))
                        break
    
    # Sort pairs: first by column (left to right), then by row (top to bottom)
    pairs.sort(key=lambda x: (x['col'], x['row']))
    
    # Then sort by row within left and right groups
    left_pairs = [p for p in pairs if p['col'] < cols//2]
    right_pairs = [p for p in pairs if p['col'] >= cols//2]
    left_pairs.sort(key=lambda x: x['row'])
    right_pairs.sort(key=lambda x: x['row'])
    
    return left_pairs + right_pairs

def solve_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Find and sort pairs
    pairs = find_and_sort_pairs(grid)
    
    # Create 4x4 output grid
    output = [[5]*4 for _ in range(4)]
    
    # Place first two pairs (left side)
    if len(pairs) >= 2:
        pair1, pair2 = pairs[0], pairs[1]
        # First pair
        output[0][0] = pair1['num']
        output[0][1] = 5 if pair1['type'] == 'horizontal' else pair1['num']
        if pair1['type'] == 'vertical':
            output[1][0] = pair1['num']
            output[1][1] = 5
        
        # Second pair
        output[2][0] = pair2['num']
        output[2][1] = 5 if pair2['type'] == 'horizontal' else pair2['num']
        if pair2['type'] == 'vertical':
            output[3][0] = pair2['num']
            output[3][1] = 5
    
    # Place second two pairs (right side)
    if len(pairs) >= 4:
        pair3, pair4 = pairs[2], pairs[3]
        # Third pair
        output[0][2] = pair3['num']
        output[0][3] = 5 if pair3['type'] == 'horizontal' else pair3['num']
        if pair3['type'] == 'vertical':
            output[1][2] = pair3['num']
            output[1][3] = 5
        
        # Fourth pair
        output[2][2] = pair4['num']
        output[2][3] = 5 if pair4['type'] == 'horizontal' else pair4['num']
        if pair4['type'] == 'vertical':
            output[3][2] = pair4['num']
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
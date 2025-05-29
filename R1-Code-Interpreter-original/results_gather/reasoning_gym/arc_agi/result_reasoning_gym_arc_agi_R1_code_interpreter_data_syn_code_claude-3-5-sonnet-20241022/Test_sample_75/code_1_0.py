def find_pairs(grid):
    pairs = []
    visited = set()
    rows, cols = len(grid), len(grid[0])
    
    # Find all pairs
    for r in range(rows):
        for c in range(cols):
            if (r, c) not in visited and grid[r][c] != 5:
                num = grid[r][c]
                # Check right and down
                for dr, dc in [(0,1), (1,0)]:
                    nr, nc = r + dr, c + dc
                    if (0 <= nr < rows and 0 <= nc < cols and 
                        grid[nr][nc] == num):
                        pairs.append({
                            'num': num,
                            'type': 'horizontal' if dr == 0 else 'vertical',
                            'row': r,
                            'col': c
                        })
                        visited.add((r, c))
                        visited.add((nr, nc))
                        break
    return pairs

def solve_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Find all pairs
    pairs = find_pairs(grid)
    
    # Sort pairs: first by row position, then by column position
    pairs.sort(key=lambda x: (x['row'], x['col']))
    
    # Create 4x4 output grid
    output = [[5]*4 for _ in range(4)]
    
    # Place pairs in output grid
    row_idx = 0
    for pair in pairs:
        # Determine placement in output grid
        if row_idx < 4:
            # Place the pair
            if pair['type'] == 'horizontal':
                output[row_idx][0] = pair['num']
                output[row_idx][1] = 5
            else:  # vertical
                output[row_idx][0] = pair['num']
                output[row_idx][1] = pair['num']
            row_idx += 1
    
    # Place the rightmost pairs
    row_idx = 0
    for pair in pairs:
        if row_idx < 4:
            if pair['type'] == 'horizontal':
                output[row_idx][2] = pair['num']
                output[row_idx][3] = 5
            else:  # vertical
                output[row_idx][2] = pair['num']
                output[row_idx][3] = pair['num']
            row_idx += 1
    
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
def find_all_pairs(grid):
    rows, cols = len(grid), len(grid[0])
    pairs = []
    visited = set()
    
    def is_valid(r, c):
        return 0 <= r < rows and 0 <= c < cols
    
    for r in range(rows):
        for c in range(cols):
            if (r, c) not in visited and grid[r][c] != 5:
                num = grid[r][c]
                # Check right and down
                for dr, dc in [(0,1), (1,0)]:
                    nr, nc = r + dr, c + dc
                    if is_valid(nr, nc) and grid[nr][nc] == num:
                        pairs.append({
                            'num': num,
                            'row': r,
                            'col': c,
                            'type': 'horizontal' if dr == 0 else 'vertical',
                            'position': 'left' if c < cols//2 else 'right'
                        })
                        visited.add((r, c))
                        visited.add((nr, nc))
                        break
    return pairs

def solve_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Find all pairs
    pairs = find_all_pairs(grid)
    
    # Sort pairs by row position
    pairs.sort(key=lambda x: (x['row'], x['col']))
    
    # Create 4x4 output grid
    output = [[5]*4 for _ in range(4)]
    
    # Process pairs
    left_pairs = [p for p in pairs if p['position'] == 'left']
    right_pairs = [p for p in pairs if p['position'] == 'right']
    
    # Place pairs in output grid
    for i, pair in enumerate(pairs):
        row = i // 2  # Determines which row to place the pair
        if i % 2 == 0:  # Left side
            output[row][0] = pair['num']
            output[row][1] = 5 if pair['type'] == 'horizontal' else pair['num']
        else:  # Right side
            output[row][2] = pair['num']
            output[row][3] = 5 if pair['type'] == 'horizontal' else pair['num']
    
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
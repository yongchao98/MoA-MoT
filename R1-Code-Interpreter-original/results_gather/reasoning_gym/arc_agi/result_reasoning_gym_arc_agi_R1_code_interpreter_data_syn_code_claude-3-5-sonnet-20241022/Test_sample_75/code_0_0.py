def find_pairs(grid):
    pairs = []
    visited = set()
    rows, cols = len(grid), len(grid[0])
    
    # Helper function to check if a position is valid
    def is_valid(r, c):
        return 0 <= r < rows and 0 <= c < cols
    
    # Find all pairs
    for r in range(rows):
        for c in range(cols):
            if (r, c) not in visited and grid[r][c] != 5:
                num = grid[r][c]
                # Check adjacent positions
                directions = [(0,1), (1,0)]  # right and down
                for dr, dc in directions:
                    nr, nc = r + dr, c + dc
                    if is_valid(nr, nc) and grid[nr][nc] == num:
                        pairs.append({
                            'num': num,
                            'pos': (r, c),
                            'second_pos': (nr, nc)
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
    
    # Sort pairs by position (top pairs first, then left pairs)
    pairs.sort(key=lambda x: (x['pos'][0], x['pos'][1]))
    
    # Create 4x4 output grid
    output = [[0]*4 for _ in range(4)]
    
    # Place pairs in output grid
    for i, pair in enumerate(pairs):
        row = i // 2  # Two pairs per row
        col = (i % 2) * 2  # Alternate between columns 0 and 2
        output[row][col] = pair['num']
        output[row][col+1] = pair['num']
        if pair['pos'][0] == pair['second_pos'][0]:  # horizontal pair
            output[row][col+1] = 5
    
    # Format output
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
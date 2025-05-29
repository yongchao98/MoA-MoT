def transform_grid(input_grid):
    n = len(input_grid)
    output = [[0]*n for _ in range(n)]
    
    def find_pattern(start_i, start_j):
        if input_grid[start_i][start_j] == 0:
            return None
            
        frame_num = input_grid[start_i][start_j]
        frame_positions = set()
        min_i = max_i = start_i
        min_j = max_j = start_j
        
        # Find pattern boundaries and frame positions
        stack = [(start_i, start_j)]
        visited = set()
        
        while stack:
            i, j = stack.pop()
            if (i, j) in visited or input_grid[i][j] == 0:
                continue
                
            visited.add((i, j))
            if input_grid[i][j] == frame_num:
                frame_positions.add((i, j))
                
            min_i = min(min_i, i)
            max_i = max(max_i, i)
            min_j = min(min_j, j)
            max_j = max(max_j, j)
            
            for ni, nj in [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]:
                if 0 <= ni < n and 0 <= nj < n and input_grid[ni][nj] != 0:
                    stack.append((ni, nj))
        
        return min_i, max_i, min_j, max_j, frame_num, frame_positions
    
    # Find and transform all patterns
    visited = set()
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] != 0 and (i,j) not in visited:
                pattern = find_pattern(i, j)
                if pattern:
                    min_i, max_i, min_j, max_j, frame_num, frame_positions = pattern
                    visited.update(frame_positions)
                    
                    # Expand boundaries by 1
                    new_min_i = max(0, min_i - 1)
                    new_max_i = min(n-1, max_i + 1)
                    new_min_j = max(0, min_j - 1)
                    new_max_j = min(n-1, max_j + 1)
                    
                    # Fill the expanded pattern
                    for r in range(new_min_i, new_max_i + 1):
                        for c in range(new_min_j, new_max_j + 1):
                            if r == new_min_i or r == new_max_i or c == new_min_j or c == new_max_j:
                                # Outer border
                                output[r][c] = frame_num
                            elif (r, c) in frame_positions:
                                # Original frame positions
                                output[r][c] = frame_num
                            else:
                                # Inner area
                                inner_num = None
                                for tr in range(min_i, max_i + 1):
                                    for tc in range(min_j, max_j + 1):
                                        if input_grid[tr][tc] != frame_num:
                                            inner_num = input_grid[tr][tc]
                                            break
                                    if inner_num:
                                        break
                                output[r][c] = inner_num if inner_num else frame_num

    return output

# Test input
test_input = [
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,6,6,6,0,0,0,0,0,0,0,0],
    [0,6,4,6,0,0,0,0,0,0,0,0],
    [0,6,6,6,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,6,6,6,6,0,0],
    [0,0,0,0,0,0,6,4,4,6,0,0],
    [0,0,0,0,0,0,6,4,4,6,0,0],
    [0,0,0,0,0,0,6,6,6,6,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0]
]

result = transform_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))
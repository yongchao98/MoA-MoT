def find_largest_block(grid):
    from collections import defaultdict, deque
    
    def bfs(start, value):
        queue = deque([start])
        visited.add(start)
        size = 0
        while queue:
            x, y = queue.popleft()
            size += 1
            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nx, ny = x + dx, y + dy
                if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited and grid[nx][ny] == value:
                    visited.add((nx, ny))
                    queue.append((nx, ny))
        return size

    visited = set()
    block_sizes = defaultdict(int)
    
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] != 3 and (i, j) not in visited:
                block_size = bfs((i, j), grid[i][j])
                block_sizes[grid[i][j]] = max(block_sizes[grid[i][j]], block_size)
    
    # Find the number with the largest block size
    largest_block_value = max(block_sizes, key=block_sizes.get)
    return largest_block_value

# Test input grid
input_grid = [
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3, 3, 3, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 2, 2, 2, 2, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3],
    [3, 6, 6, 6, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3],
    [3, 6, 6, 6, 5, 5, 5, 5, 3, 2, 2, 2, 2, 3, 3],
    [3, 6, 6, 6, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 5, 5, 5, 5, 3, 3, 4, 4, 4, 4, 3],
    [3, 3, 3, 3, 5, 5, 5, 5, 3, 3, 4, 3, 3, 4, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 3]
]

# Find the output for the test input grid
output = find_largest_block(input_grid)
print(output)
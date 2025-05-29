def find_most_significant_pattern(grid):
    def get_surrounding(i, j, val):
        pattern = []
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                ni, nj = i + di, j + dj
                if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]):
                    if grid[ni][nj] == val:
                        pattern.append((di, dj))
        return frozenset(pattern)
    
    # Count patterns for each number
    patterns = {}
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] != 7:
                val = grid[i][j]
                pattern = get_surrounding(i, j, val)
                if (val, pattern) not in patterns:
                    patterns[(val, pattern)] = 0
                patterns[(val, pattern)] += 1
    
    # Find the most common pattern
    if patterns:
        val, pattern = max(patterns.items(), key=lambda x: x[1])[0]
        return val
    return None

# Parse input
input_str = """7 7 7 7 7 7 7 7 7 9 9 7 7 7
7 7 7 7 7 7 7 7 9 9 9 7 7 7
7 7 7 7 7 7 7 7 7 9 7 7 7 7
7 7 5 7 5 7 7 7 7 7 7 7 7 7
7 7 7 5 7 7 7 7 7 7 7 7 9 9
7 7 5 5 5 7 7 7 7 7 7 9 9 9
7 7 7 7 7 7 7 7 8 7 7 7 9 7
7 7 7 7 7 7 7 8 7 8 7 7 7 7
7 7 7 7 7 7 7 7 8 7 7 7 7 7
7 7 7 9 9 7 7 7 7 7 7 5 7 5
7 7 9 9 9 7 7 7 7 7 7 7 5 7
7 7 7 9 7 7 7 7 7 7 7 5 5 5
7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 7 7 7 7 7 7 7 7 7 7 7 7"""

grid = [[int(x) for x in row.split()] for row in input_str.split('\n')]
val = find_most_significant_pattern(grid)
print(f"Most significant value: {val}")
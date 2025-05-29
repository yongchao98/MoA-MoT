def find_significant_numbers(grid):
    from collections import Counter
    import numpy as np

    # Convert grid to numpy array for easier manipulation
    grid = np.array(grid)
    rows, cols = grid.shape

    # Step 1: Frequency Analysis
    frequency = Counter(grid.flatten())

    # Step 2: Cluster Detection
    def is_cluster(r, c, num):
        # Check if the number forms a cluster by looking at adjacent cells
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols and grid[nr, nc] == num:
                return True
        return False

    clusters = Counter()
    for r in range(rows):
        for c in range(cols):
            num = grid[r, c]
            if is_cluster(r, c, num):
                clusters[num] += 1

    # Step 3: Determine Significant Numbers
    # A number is significant if it has a high frequency or forms a cluster
    significant_numbers = set()
    for num, count in frequency.items():
        if count > 1 or clusters[num] > 0:
            significant_numbers.add(num)

    # Step 4: Generate Output
    # Sort significant numbers and repeat them to form the output grid
    significant_numbers = sorted(significant_numbers)
    output_grid = [significant_numbers for _ in range(rows - 1)]

    return output_grid

# Test input grid
input_grid = [
    [1, 8, 1, 1, 1, 1, 1, 1, 1, 1],
    [8, 8, 8, 1, 1, 1, 1, 6, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 6, 6, 1],
    [1, 1, 1, 1, 1, 1, 1, 6, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 3, 3, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 9, 9, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 9, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 9, 1, 1, 1]
]

output_grid = find_significant_numbers(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
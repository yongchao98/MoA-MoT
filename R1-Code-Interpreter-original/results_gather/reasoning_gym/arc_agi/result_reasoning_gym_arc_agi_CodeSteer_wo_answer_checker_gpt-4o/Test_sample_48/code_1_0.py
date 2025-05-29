class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))
        self.rank = [1] * size
        self.size = [1] * size

    def find(self, p):
        if self.parent[p] != p:
            self.parent[p] = self.find(self.parent[p])
        return self.parent[p]

    def union(self, p, q):
        rootP = self.find(p)
        rootQ = self.find(q)
        if rootP != rootQ:
            if self.rank[rootP] > self.rank[rootQ]:
                self.parent[rootQ] = rootP
                self.size[rootP] += self.size[rootQ]
            elif self.rank[rootP] < self.rank[rootQ]:
                self.parent[rootP] = rootQ
                self.size[rootQ] += self.size[rootP]
            else:
                self.parent[rootQ] = rootP
                self.size[rootP] += self.size[rootQ]
                self.rank[rootP] += 1

    def component_size(self, p):
        rootP = self.find(p)
        return self.size[rootP]

def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    uf = UnionFind(rows * cols)

    # Directions for right and down (to avoid double counting)
    directions = [(0, 1), (1, 0)]

    # Union adjacent 1s
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 1:
                for dr, dc in directions:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols and input_grid[nr][nc] == 1:
                        uf.union(r * cols + c, nr * cols + nc)

    # Create output grid
    output_grid = [row[:] for row in input_grid]

    # Transform components of size >= 3
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 1:
                if uf.component_size(r * cols + c) >= 3:
                    output_grid[r][c] = 5

    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 6, 6, 6, 6],
    [1, 1, 1, 1, 6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1],
    [1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [1, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 1, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
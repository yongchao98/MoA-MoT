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
                self.rank[rootP] += 1
                self.size[rootP] += self.size[rootQ]

    def component_size(self, p):
        return self.size[self.find(p)]

def transform_grid(input_grid):
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    uf = UnionFind(rows * cols)

    # Directions for adjacent cells (right, down)
    directions = [(0, 1), (1, 0)]

    # Union adjacent 4s
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 4:
                for dr, dc in directions:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols and grid[nr][nc] == 4:
                        uf.union(r * cols + c, nr * cols + nc)

    # Find the largest component of 4s
    largest_component_root = None
    largest_size = 0
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 4:
                root = uf.find(r * cols + c)
                size = uf.component_size(root)
                if size > largest_size:
                    largest_size = size
                    largest_component_root = root

    # Create the output grid by replacing the largest block of 4s with 5s
    output_grid = [row[:] for row in grid]
    if largest_component_root is not None:
        for r in range(rows):
            for c in range(cols):
                if grid[r][c] == 4 and uf.find(r * cols + c) == largest_component_root:
                    output_grid[r][c] = 5

    # Print the output grid
    for row in output_grid:
        print(' '.join(map(str, row)))

# Test input grid
input_grid = """
3 3 3 3 3 3 3 3 3 3 3 5 3 3 3 3 3 3 5 3 3 3 3
3 3 3 3 4 4 4 3 3 3 3 3 3 3 4 4 4 3 3 3 3 3 3
3 3 3 3 4 4 4 3 3 3 3 3 3 3 4 4 4 3 4 4 4 4 3
3 3 3 3 4 4 4 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3
3 3 3 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 4 4 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3
3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3
3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4
3 3 4 4 4 4 4 3 3 3 3 4 4 4 3 3 3 3 3 4 4 4 4
3 3 4 4 4 4 4 3 3 3 3 4 4 4 3 3 3 3 3 4 4 4 4
3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 4 4 4 3 3 3
5 3 4 4 4 4 4 3 3 4 4 4 4 3 3 3 3 4 4 4 3 3 5
3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 4 4 4 3 3 3
3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 5 3 3 3 3 3 3 5 3 3 3 3
"""

transform_grid(input_grid)
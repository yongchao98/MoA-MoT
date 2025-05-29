class Node:
    def __init__(self):
        self.left = self.right = self.up = self.down = self
        self.column = None

class ColumnNode(Node):
    def __init__(self, name):
        super().__init__()
        self.size = 0
        self.name = name

class DancingLinks:
    def __init__(self, matrix):
        self.header = ColumnNode("header")
        self.columns = []
        self.solution = []

        # Create column nodes
        for i in range(len(matrix[0])):
            column = ColumnNode(i)
            self.columns.append(column)
            self.header = self._link_horizontal(self.header, column)

        # Create data nodes
        for i, row in enumerate(matrix):
            prev = None
            for j, cell in enumerate(row):
                if cell == 1:
                    column = self.columns[j]
                    node = Node()
                    node.column = column
                    column.size += 1
                    column = self._link_vertical(column, node)
                    if prev is None:
                        prev = node
                    else:
                        prev = self._link_horizontal(prev, node)

    def _link_horizontal(self, left, right):
        left.right = right
        right.left = left
        return right

    def _link_vertical(self, up, down):
        up.down = down
        down.up = up
        return down

    def cover(self, column):
        column.right.left = column.left
        column.left.right = column.right
        for row in self._iterate(column.down, 'down'):
            for node in self._iterate(row.right, 'right'):
                node.down.up = node.up
                node.up.down = node
                node.column.size -= 1

    def uncover(self, column):
        for row in self._iterate(column.up, 'up'):
            for node in self._iterate(row.left, 'left'):
                node.column.size += 1
                node.down.up = node
                node.up.down = node
        column.right.left = column
        column.left.right = column

    def _iterate(self, start, direction):
        node = start
        while True:
            yield node
            node = getattr(node, direction)
            if node == start:
                break

    def search(self, k=0):
        if self.header.right == self.header:
            return True

        column = min(self._iterate(self.header.right, 'right'), key=lambda c: c.size)
        self.cover(column)

        for row in self._iterate(column.down, 'down'):
            self.solution.append(row)
            for node in self._iterate(row.right, 'right'):
                self.cover(node.column)

            if self.search(k + 1):
                return True

            self.solution.pop()
            for node in self._iterate(row.left, 'left'):
                self.uncover(node.column)

        self.uncover(column)
        return False

def create_matrix(n, blocked_positions):
    matrix = []
    for i in range(n):
        for j in range(n):
            if (i, j) in blocked_positions:
                continue
            row = [0] * (6 * n - 2)
            row[i] = 1
            row[n + j] = 1
            row[2 * n + i + j] = 1
            row[5 * n - 2 + i - j] = 1
            matrix.append(row)
    return matrix

def solve_n_queens_dancing_links(n, blocked_positions):
    matrix = create_matrix(n, blocked_positions)
    print("Matrix constructed with dimensions:", len(matrix), "x", len(matrix[0]) if matrix else 0)
    for row in matrix:
        print(row)
    dlx = DancingLinks(matrix)
    if dlx.search():
        result = []
        for row in dlx.solution:
            for node in dlx._iterate(row.right, 'right'):
                if node.column.name < n:
                    result.append((node.column.name, node.column.name))
        print("<<<" + ", ".join(f"{r} {c}" for r, c in result) + ">>>")
    else:
        print("No solution found")

blocked_positions = {(1, 6), (2, 0)}
solve_n_queens_dancing_links(8, blocked_positions)
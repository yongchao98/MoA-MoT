class Node:
    def __init__(self, row=None, col=None):
        self.left = self.right = self.up = self.down = self
        self.row = row
        self.col = col

class ColumnNode(Node):
    def __init__(self, name):
        super().__init__()
        self.size = 0
        self.name = name

class DLX:
    def __init__(self, matrix):
        self.header = ColumnNode("header")
        self.columns = []
        self.solution = []
        self.build_linked_list(matrix)

    def build_linked_list(self, matrix):
        num_cols = len(matrix[0])
        column_nodes = [ColumnNode(i) for i in range(num_cols)]
        self.columns = column_nodes

        # Link header to column nodes
        current = self.header
        for col_node in column_nodes:
            current.right = col_node
            col_node.left = current
            current = col_node
        current.right = self.header
        self.header.left = current

        # Create nodes for each 1 in the matrix
        for r, row in enumerate(matrix):
            prev_node = None
            for c, val in enumerate(row):
                if val == 1:
                    col_node = column_nodes[c]
                    new_node = Node(r, c)
                    col_node.size += 1

                    # Link vertically
                    new_node.down = col_node
                    new_node.up = col_node.up
                    col_node.up.down = new_node
                    col_node.up = new_node

                    # Link horizontally
                    if prev_node is None:
                        prev_node = new_node
                    else:
                        prev_node.right = new_node
                        new_node.left = prev_node
                        prev_node = new_node

            if prev_node:
                first_node = prev_node.right
                prev_node.right = first_node
                first_node.left = prev_node

    def cover(self, col_node):
        col_node.right.left = col_node.left
        col_node.left.right = col_node.right
        for row_node in self.iterate(col_node.down, 'down'):
            for node in self.iterate(row_node.right, 'right'):
                node.down.up = node.up
                node.up.down = node.down
                self.columns[node.col].size -= 1

    def uncover(self, col_node):
        for row_node in self.iterate(col_node.up, 'up'):
            for node in self.iterate(row_node.left, 'left'):
                self.columns[node.col].size += 1
                node.down.up = node
                node.up.down = node
        col_node.right.left = col_node
        col_node.left.right = col_node

    def iterate(self, start, direction):
        node = start
        while node != start or node == start:
            yield node
            node = getattr(node, direction)

    def search(self, k=0):
        if self.header.right == self.header:
            return True

        # Choose the column with the smallest size
        col_node = min(self.iterate(self.header.right, 'right'), key=lambda c: c.size)
        self.cover(col_node)

        for row_node in self.iterate(col_node.down, 'down'):
            self.solution.append(row_node.row)
            for node in self.iterate(row_node.right, 'right'):
                self.cover(self.columns[node.col])

            if self.search(k + 1):
                return True

            self.solution.pop()
            for node in self.iterate(row_node.left, 'left'):
                self.uncover(self.columns[node.col])

        self.uncover(col_node)
        return False

def solve_puzzle_with_dlx():
    # Define the constraints and the matrix for the DLX algorithm
    # This part involves creating the exact cover matrix for the problem
    # For simplicity, this example assumes the matrix is already defined
    # You would need to create this matrix based on the problem constraints

    # Example matrix (this needs to be constructed based on the problem)
    matrix = [
        # This is a placeholder; the actual matrix needs to be constructed
    ]

    dlx = DLX(matrix)
    if dlx.search():
        # Convert the solution back to the grid format
        # This part involves interpreting the solution from the DLX algorithm
        # and printing it in the required format
        pass

solve_puzzle_with_dlx()
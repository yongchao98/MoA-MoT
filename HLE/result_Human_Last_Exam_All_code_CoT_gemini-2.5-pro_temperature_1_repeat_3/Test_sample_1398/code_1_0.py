import sys

def solve():
    """
    This function solves the grid transformation puzzle by learning the rules from examples
    and applying them to a test grid.
    The rule is assumed to be a function of the cell's value and the sums of its
    Moore and Von Neumann neighbors, with wrap-around boundaries.
    """
    
    # Input Grids from the problem description
    inputs = [
        [[0, 1, 0, 0, 1], [1, 0, 0, 1, 0], [0, 0, 1, 0, 0], [1, 0, 0, 0, 1], [0, 1, 0, 1, 0]],
        [[1, 1, 0, 1, 0], [0, 0, 1, 0, 1], [1, 0, 0, 1, 0], [0, 1, 1, 0, 0], [1, 0, 0, 1, 1]],
        [[0, 0, 1, 1, 0], [1, 0, 0, 0, 1], [0, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 1, 0, 0, 1]],
        [[1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1]],
        [[0, 0, 0, 0, 0], [0, 1, 1, 1, 0], [0, 1, 0, 1, 0], [0, 1, 1, 1, 0], [0, 0, 0, 0, 0]],
    ]
    
    # Output Grids from the problem description
    outputs = [
        [[1, 0, 1, 1, 0], [0, 1, 1, 1, 1], [1, 1, 0, 1, 1], [0, 1, 1, 1, 0], [1, 0, 1, 0, 1]],
        [[0, 1, 1, 1, 1], [1, 0, 1, 0, 1], [0, 0, 0, 1, 1], [1, 1, 1, 0, 1], [0, 1, 1, 1, 0]],
        [[0, 1, 0, 1, 1], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], [1, 0, 0, 1, 1], [1, 0, 1, 1, 0]],
        [[0, 1, 1, 1, 0], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [0, 1, 1, 1, 0]],
        [[0, 1, 1, 1, 0], [1, 1, 0, 1, 1], [1, 0, 0, 0, 1], [1, 1, 0, 1, 1], [0, 1, 1, 1, 0]],
    ]
    
    # The 6th test input grid
    test_input = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    def get_neighbor_sums(grid, r, c):
        """Calculates the Moore and Von Neumann neighbor sums for a cell."""
        rows, cols = len(grid), len(grid[0])
        moore_sum = 0
        vn_sum = 0
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                # Apply wrap-around (toroidal) boundary conditions
                nr, nc = (r + dr + rows) % rows, (c + dc + cols) % cols
                moore_sum += grid[nr][nc]
                # Von Neumann neighbors are up, down, left, right
                if abs(dr) + abs(dc) == 1:
                    vn_sum += grid[nr][nc]
        return moore_sum, vn_sum

    # This dictionary will store the learned rules
    rules = {}
    
    # Learn the rules from the 5 examples
    for i in range(len(inputs)):
        input_grid = inputs[i]
        output_grid = outputs[i]
        rows, cols = len(input_grid), len(input_grid[0])
        for r in range(rows):
            for c in range(cols):
                val = input_grid[r][c]
                m_sum, v_sum = get_neighbor_sums(input_grid, r, c)
                key = (val, m_sum, v_sum)
                out = output_grid[r][c]
                
                # Check for contradictions
                if key in rules and rules[key] != out:
                    print(f"Error: Contradiction found for key {key}. Rule is not consistent.")
                    sys.exit(1)
                rules[key] = out

    # Apply the learned rules to the test grid
    test_output = [[0 for _ in range(5)] for _ in range(5)]
    rows, cols = len(test_input), len(test_input[0])
    for r in range(rows):
        for c in range(cols):
            val = test_input[r][c]
            m_sum, v_sum = get_neighbor_sums(test_input, r, c)
            key = (val, m_sum, v_sum)
            
            if key in rules:
                test_output[r][c] = rules[key]
            else:
                # This case means the test grid has a configuration not seen in the examples.
                print(f"Error: No rule found for key {key}. Cannot determine output.")
                # As a fallback, you might assume a default output, e.g., the input value itself.
                # However, for this problem, it's expected that all rules can be learned.
                test_output[r][c] = -1 # Indicate an error

    # Print the resulting output grid
    final_string = ""
    for r in range(rows):
        row_str = " ".join(map(str, test_output[r]))
        final_string += "".join(map(str, test_output[r]))
        print(row_str)

# Execute the solver
solve()
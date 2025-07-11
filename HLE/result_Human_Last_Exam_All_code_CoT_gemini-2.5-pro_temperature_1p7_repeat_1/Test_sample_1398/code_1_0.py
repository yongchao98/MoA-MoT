import numpy as np

def solve_grid_puzzle():
    """
    Solves the grid puzzle by learning a cellular automaton rule from examples
    and applying it to the test grid.
    """
    
    # Define the 5 example grids and the 1 test grid
    grids = {
        '1': {
            'input': np.array([
                [0, 1, 0, 0, 1], [1, 0, 0, 1, 0], [0, 0, 1, 0, 0], [1, 0, 0, 0, 1], [0, 1, 0, 1, 0]
            ]),
            'output': np.array([
                [1, 0, 1, 1, 0], [0, 1, 1, 1, 1], [1, 1, 0, 1, 1], [0, 1, 1, 1, 0], [1, 0, 1, 0, 1]
            ])
        },
        '2': {
            'input': np.array([
                [1, 1, 0, 1, 0], [0, 0, 1, 0, 1], [1, 0, 0, 1, 0], [0, 1, 1, 0, 0], [1, 0, 0, 1, 1]
            ]),
            'output': np.array([
                [0, 1, 1, 1, 1], [1, 0, 1, 0, 1], [0, 0, 0, 1, 1], [1, 1, 1, 0, 1], [0, 1, 1, 1, 0]
            ])
        },
        '3': {
            'input': np.array([
                [0, 0, 1, 1, 0], [1, 0, 0, 0, 1], [0, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 1, 0, 0, 1]
            ]),
            'output': np.array([
                [0, 1, 0, 1, 1], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], [1, 0, 0, 1, 1], [1, 0, 1, 1, 0]
            ])
        },
        '4': {
            'input': np.array([
                [1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1]
            ]),
            'output': np.array([
                [0, 1, 1, 1, 0], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [0, 1, 1, 1, 0]
            ])
        },
        '5': {
            'input': np.array([
                [0, 0, 0, 0, 0], [0, 1, 1, 1, 0], [0, 1, 0, 1, 0], [0, 1, 1, 1, 0], [0, 0, 0, 0, 0]
            ]),
            'output': np.array([
                [0, 1, 1, 1, 0], [1, 1, 0, 1, 1], [1, 0, 0, 0, 1], [1, 1, 0, 1, 1], [0, 1, 1, 1, 0]
            ])
        }
    }

    test_input = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ])

    # This dictionary will store the learned rules: (input_val, neighbor_sum) -> output_val
    rules = {}

    def get_neighbor_sum(grid, r, c):
        """Calculates the sum of the 8 neighbors for a cell (r, c)."""
        s = 0
        rows, cols = grid.shape
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    s += grid[nr, nc]
        return s

    # Learn the rules from the 5 examples
    for i in range(1, 6):
        input_grid = grids[str(i)]['input']
        output_grid = grids[str(i)]['output']
        for r in range(5):
            for c in range(5):
                val = input_grid[r, c]
                n_sum = get_neighbor_sum(input_grid, r, c)
                out = output_grid[r, c]
                
                key = (val, n_sum)
                # Store the rule, checking for contradictions
                if key in rules and rules[key] != out:
                    # This check confirms the hypothesis that the rule is consistent
                    raise Exception(f"Contradiction found for rule {key} -> {rules[key]} vs {out}")
                rules[key] = out
    
    # Apply the learned rules to the test grid
    final_output_grid = np.zeros((5, 5), dtype=int)
    for r in range(5):
        for c in range(5):
            val = test_input[r, c]
            n_sum = get_neighbor_sum(test_input, r, c)
            key = (val, n_sum)
            if key in rules:
                final_output_grid[r, c] = rules[key]
            else:
                # Handle cases where a rule was not learned from the training data.
                # In a robust solution, one might define a default behavior.
                # For this puzzle, all necessary rules are learnable.
                raise Exception(f"Rule not found for key {key}")
                
    # Flatten the 5x5 grid to a 25-character string and print
    result_string = "".join(map(str, final_output_grid.flatten()))
    print(result_string)

solve_grid_puzzle()
import numpy as np

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by simulating all possible totalistic rules
    and matching the generated patterns to the given states.
    """
    # Step 1: Define the 9 states as numpy arrays
    states = {}
    states[1] = np.array([
        [0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]
    ])
    states[2] = np.array([
        [1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]
    ])
    states[3] = np.array([
        [0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]
    ])
    states[4] = np.array([
        [0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]
    ])
    states[5] = np.array([
        [1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]
    ])
    states[6] = np.array([
        [1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]
ax])
    states[7] = np.array([
        [1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]
    ])
    states[8] = np.array([
        [0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]
    ])
    states[9] = np.array([
        [1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]
    ])

    # Step 2: Helper functions
    def trim_grid(grid):
        if np.sum(grid) == 0: return np.array([[0]])
        rows, cols = np.any(grid, axis=1), np.any(grid, axis=0)
        if not np.any(rows) or not np.any(cols): return np.array([[0]])
        rmin, rmax = np.where(rows)[0][[0, -1]]
        cmin, cmax = np.where(cols)[0][[0, -1]]
        return grid[rmin:rmax+1, cmin:cmax+1]

    trimmed_states = {i: trim_grid(s) for i, s in states.items()}

    def evolve(grid, rule_set):
        padded_grid = np.pad(grid, 1, 'constant', constant_values=0)
        new_grid = np.zeros_like(padded_grid)
        for r in range(1, padded_grid.shape[0] - 1):
            for c in range(1, padded_grid.shape[1] - 1):
                neighborhood_sum = np.sum(padded_grid[r-1:r+2, c-1:c+2])
                if neighborhood_sum in rule_set:
                    new_grid[r, c] = 1
        return new_grid

    def find_match(gen_grid, available_labels):
        trimmed_gen = trim_grid(gen_grid)
        for label in available_labels:
            if trimmed_gen.shape == trimmed_states[label].shape and \
               np.array_equal(trimmed_gen, trimmed_states[label]):
                return label
        return None

    # Step 3 & 4: Main search and match loop
    found_rules = {}
    found_labels = set()
    all_labels = set(range(1, 10))
    
    for r_idx in range(1024):
        rule_set = {i for i, bit in enumerate(bin(r_idx)[2:].zfill(10)[::-1]) if bit == '1'}
        if not rule_set: continue

        grids = [np.array([[1]])]
        for _ in range(6):
            grids.append(evolve(grids[-1], rule_set))

        time_steps_map = {1: (2, 3, 4), 2: (3, 4, 5), 3: (4, 5, 6)}

        for rule_num, (t1, t2, t3) in time_steps_map.items():
            if rule_num in found_rules: continue

            current_available_labels = all_labels - found_labels
            
            l1 = find_match(grids[t1], current_available_labels)
            if l1 is None: continue
            
            current_available_labels.remove(l1)
            l2 = find_match(grids[t2], current_available_labels)
            if l2 is None: continue
            
            current_available_labels.remove(l2)
            l3 = find_match(grids[t3], current_available_labels)
            if l3 is None: continue

            found_rules[rule_num] = (l1, l2, l3)
            found_labels.update([l1, l2, l3])
            break
            
        if len(found_rules) == 3:
            break

    # Step 5: Construct and print the final answer
    r1_labels = "".join(map(str, found_rules[1]))
    r2_labels = "".join(map(str, found_rules[2]))
    r3_labels = "".join(map(str, found_rules[3]))
    
    result_string = f"{{{r1_labels},{r2_labels},{r3_labels}}}"
    print(result_string)

solve_ca_puzzle()
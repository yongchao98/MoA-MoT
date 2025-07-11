import numpy as np
from itertools import permutations

def parse_state(s):
    return np.array([list(map(int, row.split())) for row in s.strip().split('\n')])

# --- State Definitions ---
states_str = {
    1: """
0 0 0 1 1 1 1 1 0 0 0
0 0 1 0 0 1 0 0 1 0 0
0 1 1 1 1 0 1 1 1 1 0
1 0 1 0 0 1 0 0 1 0 1
1 0 1 0 0 1 0 0 1 0 1
1 1 0 1 1 0 1 1 0 1 1
1 0 1 0 0 1 0 0 1 0 1
1 0 1 0 0 1 0 0 1 0 1
0 1 1 1 1 0 1 1 1 1 0
0 0 1 0 0 1 0 0 1 0 0
0 0 0 1 1 1 1 1 0 0 0
""",
    2: """
1 0 0 0 1 0 0 0 1 0 0 0 1
0 1 1 1 0 0 0 0 0 1 1 1 0
0 1 1 1 0 0 0 0 0 1 1 1 0
0 1 1 1 1 1 1 1 1 1 1 1 0
1 0 0 1 1 1 1 1 1 1 0 0 1
0 0 0 1 1 1 1 1 1 1 0 0 0
0 0 0 1 1 1 1 1 1 1 0 0 0
0 0 0 1 1 1 1 1 1 1 0 0 0
1 0 0 1 1 1 1 1 1 1 0 0 1
0 1 1 1 1 1 1 1 1 1 1 1 0
0 1 1 1 0 0 0 0 0 1 1 1 0
0 1 1 1 0 0 0 0 0 1 1 1 0
1 0 0 0 1 0 0 0 1 0 0 0 1
""",
    3: """
0 0 1 0 0
0 1 0 1 0
1 0 1 0 1
0 1 0 1 0
0 0 1 0 0
""",
    4: """
0 0 0 1 0 1 0 0 0
0 1 1 0 0 0 1 1 0
0 1 1 0 0 0 1 1 0
1 0 0 0 0 0 0 0 1
0 0 0 0 1 0 0 0 0
1 0 0 0 0 0 0 0 1
0 1 1 0 0 0 1 1 0
0 1 1 0 0 0 1 1 0
0 0 0 1 0 1 0 0 0
""",
    5: """
1 1 1 0 0 1 0 0 1 1 1
1 0 1 0 1 1 1 0 1 0 1
1 1 1 0 0 1 0 0 1 1 1
0 0 0 1 1 1 1 1 0 0 0
0 1 0 1 1 0 1 1 0 1 0
1 1 1 1 0 0 0 1 1 1 1
0 1 0 1 1 0 1 1 0 1 0
0 0 0 1 1 1 1 1 0 0 0
1 1 1 0 0 1 0 0 1 1 1
1 0 1 0 1 1 1 0 1 0 1
1 1 1 0 0 1 0 0 1 1 1
""",
    6: """
1 0 0 0 0 0 0 0 1
0 1 1 0 0 0 1 1 0
0 1 1 1 0 1 1 1 0
0 0 1 1 1 1 1 0 0
0 0 0 1 1 1 0 0 0
0 0 1 1 1 1 1 0 0
0 1 1 1 0 1 1 1 0
0 1 1 0 0 0 1 1 0
1 0 0 0 0 0 0 0 1
""",
    7: """
1 1 0 0 0 1 1
1 0 0 0 0 0 1
0 0 1 1 1 0 0
0 0 1 1 1 0 0
0 0 1 1 1 0 0
1 0 0 0 0 0 1
1 1 0 0 0 1 1
""",
    8: """
0 0 0 0 1 0 0 0 0
0 0 1 1 0 1 1 0 0
0 1 1 0 0 0 1 1 0
0 1 0 1 1 1 0 1 0
1 0 0 1 0 1 0 0 1
0 1 0 1 1 1 0 1 0
0 1 1 0 0 0 1 1 0
0 0 1 1 0 1 1 0 0
0 0 0 0 1 0 0 0 0
""",
    9: """
1 0 1 0 1 0 1
0 1 0 0 0 1 0
1 0 0 1 0 0 1
0 0 1 0 1 0 0
1 0 0 1 0 0 1
0 1 0 0 0 1 0
1 0 1 0 1 0 1
"""
}

states = {k: parse_state(v) for k, v in states_str.items()}

# --- CA Helper Functions ---

def trim_grid(grid):
    if grid.sum() == 0: return np.array([[0]])
    rows = np.any(grid, axis=1)
    cols = np.any(grid, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    return grid[rmin:rmax+1, cmin:cmax+1]

def deduce_rule(grid_t0, grid_t1):
    rule = {}
    h0, w0 = grid_t0.shape
    h1, w1 = grid_t1.shape
    
    # Pad grid_t0 to calculate sums for the larger grid_t1
    pad_h = (h1 - h0) // 2 + 1
    pad_w = (w1 - w0) // 2 + 1
    padded_t0 = np.pad(grid_t0, ((pad_h, pad_h), (pad_w, pad_w)), 'constant')
    
    h_padded, w_padded = padded_t0.shape
    
    for r in range(h1):
        for c in range(w1):
            # The window in padded_t0 corresponds to the cell (r,c) in grid_t1
            r_start = r + pad_h - (h1 // 2)
            c_start = c + pad_w - (w1 // 2)
            
            s = np.sum(padded_t0[r_start-1:r_start+2, c_start-1:c_start+2])
            next_state = grid_t1[r, c]
            
            if s in rule and rule[s] != next_state:
                return None  # Conflict
            rule[s] = next_state
    return rule

def evolve(grid, rule_map):
    rule = [rule_map.get(i, 0) for i in range(10)]
    padded_grid = np.pad(grid, ((1, 1), (1, 1)), 'constant')
    h, w = padded_grid.shape
    new_grid = np.zeros_like(padded_grid)
    for r in range(h - 2):
        for c in range(w - 2):
            s = np.sum(padded_grid[r:r+3, c:c+3])
            if s < 10:
                new_grid[r+1, c+1] = rule[s]
    return trim_grid(new_grid)

def check_chain(s1_id, s2_id, s3_id):
    s1, s2, s3 = states[s1_id], states[s2_id], states[s3_id]
    
    rule12 = deduce_rule(s1, s2)
    if rule12 is None: return False
    
    rule23 = deduce_rule(s2, s3)
    if rule23 is None: return False
    
    # Merge rules
    merged_rule = rule12.copy()
    for s, state in rule23.items():
        if s in merged_rule and merged_rule[s] != state:
            return False
        merged_rule[s] = state
        
    # Verify evolution with the merged rule
    s2_gen = evolve(s1, merged_rule)
    if not np.array_equal(s2_gen, s2): return False
    
    s3_gen = evolve(s2, merged_rule)
    if not np.array_equal(s3_gen, s3): return False
    
    return True

# --- Main Logic ---

# 1. Map states to time steps by size
time_map = {
    2: [k for k, v in states.items() if v.shape == (5, 5)],
    3: [k for k, v in states.items() if v.shape == (7, 7)],
    4: [k for k, v in states.items() if v.shape == (9, 9)],
    5: [k for k, v in states.items() if v.shape == (11, 11)],
    6: [k for k, v in states.items() if v.shape == (13, 13)],
}

# 2. Use symmetry to find Rule 1
# A D4-symmetric initial state under a totalistic rule yields D4-symmetric descendants.
# States with D4 symmetry: #3, #9, #8, #2
rule1_s2 = time_map[2][0] # #3 is the only t=2 state
rule1_s3 = 9 # #9 is the only D4-symmetric t=3 state
rule1_s4 = 8 # #8 is the only D4-symmetric t=4 state
rule1_labels = f"{rule1_s2}{rule1_s3}{rule1_s4}"

# 3. Find Rules 2 and 3 by checking all remaining permutations
rem_s3 = [s for s in time_map[3] if s not in [rule1_s3]]
rem_s4 = [s for s in time_map[4] if s not in [rule1_s4]]
rem_s5 = [s for s in time_map[5]]
rem_s6 = [s for s in time_map[6]]

found_solution = False
for p3 in permutations(rem_s3):
    if found_solution: break
    for p4 in permutations(rem_s4):
        if found_solution: break
        for p5 in permutations(rem_s5):
            # Assign states to chains based on permutations
            rule2_chain = (p3[0], p4[0], p5[0])
            rule3_chain = (p4[1], p5[1], rem_s6[0])

            # Check if both chains are valid
            if check_chain(*rule2_chain) and check_chain(*rule3_chain):
                rule2_labels = "".join(map(str, rule2_chain))
                rule3_labels = "".join(map(str, rule3_chain))
                found_solution = True
                break

# 4. Print the final result
print(f"{{{rule1_labels},{rule2_labels},{rule3_labels}}}")
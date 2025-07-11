import numpy as np
import itertools
from scipy.signal import convolve2d

def solve():
    """
    Solves the cellular automaton puzzle by systematically searching for rules
    that generate the given state sequences.
    """
    states_raw = {
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

    states = {i: np.array([list(map(int, line.split())) for line in s.strip().split('\n')]) for i, s in states_raw.items()}
    
    def trim(grid):
        if grid.sum() == 0: return np.array([[0]])
        rows = np.any(grid, axis=1)
        cols = np.any(grid, axis=0)
        if not np.any(rows): return np.array([[0]])
        rmin, rmax = np.where(rows)[0][[0, -1]]
        cmin, cmax = np.where(cols)[0][[0, -1]]
        return grid[rmin:rmax+1, cmin:cmax+1]

    def compare_grids(grid1, grid2):
        trimmed1 = trim(grid1)
        trimmed2 = trim(grid2)
        if trimmed1.shape != trimmed2.shape: return False
        return np.array_equal(trimmed1, trimmed2)

    kernel = np.ones((3, 3), dtype=int)
    def evolve(grid, rule_set):
        sums = convolve2d(grid, kernel, mode='full')
        next_grid = np.zeros_like(sums, dtype=int)
        next_grid[np.isin(sums, list(rule_set))] = 1
        return next_grid

    def test_rule_sequence(rule, t_start, state_labels):
        if 1 not in rule: return False
        
        current_state = np.array([[1]])
        for _ in range(t_start):
            current_state = evolve(current_state, rule)
        
        if not compare_grids(current_state, states[state_labels[0]]): return False
        
        for i in range(len(state_labels) - 1):
            current_state = evolve(states[state_labels[i]], rule)
            if not compare_grids(current_state, states[state_labels[i+1]]): return False
            
        return True

    times_to_states = {t: [k for k, v in states.items() if v.shape[0] == 2*t + 1] for t in range(2, 7)}
    
    rule_defs = [
        (1, (2, 3, 4)),
        (2, (3, 4, 5)),
        (3, (4, 5, 6)),
    ]

    all_possible_rules = [set(r) for i in range(1, 1 << 10) for r in [itertools.compress(range(10), (i >> j) & 1 for j in range(10))]]

    solutions = {}
    found_rules = set()
    found_labels = set()

    for rule_id, time_steps in rule_defs:
        t_start, t_mid, t_end = time_steps
        candidate_sequences = itertools.product(times_to_states[t_start], times_to_states[t_mid], times_to_states[t_end])
        
        for seq in candidate_sequences:
            if any(s in found_labels for s in seq): continue
            
            for rule in all_possible_rules:
                rule_tuple = tuple(sorted(list(rule)))
                if rule_tuple in found_rules: continue

                if test_rule_sequence(rule, t_start, seq):
                    solutions[rule_id] = "".join(map(str, seq))
                    found_rules.add(rule_tuple)
                    found_labels.update(seq)
                    break
            if rule_id in solutions:
                break
    
    result_string = f"{{{solutions.get(1, '')},{solutions.get(2, '')},{solutions.get(3, '')}}}"
    print(result_string)

solve()
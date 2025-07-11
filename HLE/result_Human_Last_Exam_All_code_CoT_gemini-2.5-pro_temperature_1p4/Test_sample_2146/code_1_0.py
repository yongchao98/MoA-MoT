import numpy as np

def solve_cellular_automaton():
    """
    Solves the cellular automaton puzzle by identifying the state sequences for three distinct rules.
    """

    # Step 1: Parse the states from the problem description into numpy arrays.
    states = {
        1: np.array([[0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]]),
        2: np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]]),
        3: np.array([[0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]]),
        4: np.array([[0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]]),
        5: np.array([[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]]),
        6: np.array([[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]]),
        7: np.array([[1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]]),
        8: np.array([[0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]]),
        9: np.array([[1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]]),
    }

    # Group state labels by time step, inferred from grid size (size = 2*t + 1).
    states_by_size = {}
    for label, grid in states.items():
        size = grid.shape[0]
        if size not in states_by_size: states_by_size[size] = []
        states_by_size[size].append(label)
    
    states_by_t = {int((size - 1) / 2): labels for size, labels in states_by_size.items()}

    # --- Helper Functions ---
    def infer_rule_from_transition(grid_t, grid_t_plus_1):
        """Infers the CA rule from a transition. Returns the rule {sum: next_state} or None if inconsistent."""
        if grid_t_plus_1.shape[0] != grid_t.shape[0] + 2: return None
        
        padded_grid = np.pad(grid_t, pad_width=2, mode='constant', constant_values=0)
        rule = {}
        rows_next, cols_next = grid_t_plus_1.shape
        for r in range(rows_next):
            for c in range(cols_next):
                window = padded_grid[r:r+3, c:c+3]
                sum_val = np.sum(window)
                output_val = grid_t_plus_1[r, c]
                if sum_val in rule and rule[sum_val] != output_val: return None
                rule[sum_val] = output_val
        return rule

    def merge_rules(rule1, rule2):
        """Merges two rules, checking for consistency. Returns merged rule or None."""
        merged = rule1.copy()
        for key, value in rule2.items():
            if key in merged and merged[key] != value: return None
            merged[key] = value
        return merged

    # --- Search Logic ---
    solution = {}
    found_labels = set()

    # Find Rule 1 (t=2,3,4)
    for s2 in states_by_t[2]:
        for s3 in states_by_t[3]:
            rule1_p1 = infer_rule_from_transition(states[s2], states[s3])
            if not rule1_p1: continue
            for s4 in states_by_t[4]:
                rule1_p2 = infer_rule_from_transition(states[s3], states[s4])
                if rule1_p2 and merge_rules(rule1_p1, rule1_p2):
                    solution[1] = (s2, s3, s4)
                    found_labels.update(solution[1])
                    break
            if 1 in solution: break
        if 1 in solution: break

    # Find Rule 3 (t=4,5,6)
    t4_cands = [l for l in states_by_t[4] if l not in found_labels]
    t5_cands = [l for l in states_by_t[5] if l not in found_labels]
    t6_cands = [l for l in states_by_t[6] if l not in found_labels]
    for s4 in t4_cands:
        for s5 in t5_cands:
            rule3_p1 = infer_rule_from_transition(states[s4], states[s5])
            if not rule3_p1: continue
            for s6 in t6_cands:
                rule3_p2 = infer_rule_from_transition(states[s5], states[s6])
                if rule3_p2 and merge_rules(rule3_p1, rule3_p2):
                    solution[3] = (s4, s5, s6)
                    found_labels.update(solution[3])
                    break
            if 3 in solution: break
        if 3 in solution: break

    # Rule 2 are the remaining labels
    remaining = set(states.keys()) - found_labels
    s3 = [l for l in remaining if l in states_by_t[3]][0]
    s4 = [l for l in remaining if l in states_by_t[4]][0]
    s5 = [l for l in remaining if l in states_by_t[5]][0]
    solution[2] = (s3, s4, s5)

    # --- Format and Print Output ---
    r1_labels = "".join(map(str, solution[1]))
    r2_labels = "".join(map(str, solution[2]))
    r3_labels = "".join(map(str, solution[3]))
    
    print(f"{{{r1_labels},{r2_labels},{r3_labels}}}")

solve_cellular_automaton()
import numpy as np

def solve_ca_puzzle():
    """
    Solves the cellular automaton state matching puzzle.
    """
    # 1. Parse and represent the states as numpy arrays
    s1 = np.array([[0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]])
    s2 = np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]])
    s3 = np.array([[0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]])
    s4 = np.array([[0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]])
    s5 = np.array([[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]])
    s6 = np.array([[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]])
    s7 = np.array([[1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]])
    s8 = np.array([[0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]])
    s9 = np.array([[1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]])

    states = {1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7, 8: s8, 9: s9}
    
    # Pad states for easier processing
    padded_states = {}
    GRID_SIZE = 25
    for i, s in states.items():
        h, w = s.shape
        padded_s = np.zeros((GRID_SIZE, GRID_SIZE), dtype=int)
        start_h = (GRID_SIZE - h) // 2
        start_w = (GRID_SIZE - w) // 2
        padded_s[start_h:start_h+h, start_w:start_w+w] = s
        padded_states[i] = padded_s

    # 2. Group states by size/timestep
    size_map = {
        2: [3], 3: [7, 9], 4: [4, 6, 8], 5: [1, 5], 6: [2]
    }

    # 3. & 4. Function to check transitions and deduce rules
    def check_transition(s1_idx, s2_idx, base_rule=None):
        grid1 = padded_states[s1_idx]
        grid2 = padded_states[s2_idx]
        rule = {} if base_rule is None else base_rule.copy()
        
        # Define checking area based on pattern sizes
        check_dim = max(states[s1_idx].shape[0], states[s2_idx].shape[0]) + 2
        start = (GRID_SIZE - check_dim) // 2
        end = start + check_dim

        for r in range(start, end):
            for c in range(start, end):
                s = np.sum(grid1[r-1:r+2, c-1:c+2])
                next_state = grid2[r, c]
                if s in rule:
                    if rule[s] != next_state:
                        return None
                else:
                    rule[s] = next_state
        return rule

    # 5. Search for valid sequences for each rule
    possible_rule1s, possible_rule2s, possible_rule3s = [], [], []

    # Rule 1: t=2,3,4
    for s_t3 in size_map[3]:
        rule1_cand = check_transition(size_map[2][0], s_t3)
        if rule1_cand:
            for s_t4 in size_map[4]:
                final_rule = check_transition(s_t3, s_t4, rule1_cand)
                if final_rule:
                    possible_rule1s.append(((size_map[2][0], s_t3, s_t4), final_rule))

    # Rule 2: t=3,4,5
    for s_t3 in size_map[3]:
        for s_t4 in size_map[4]:
            rule2_cand = check_transition(s_t3, s_t4)
            if rule2_cand:
                for s_t5 in size_map[5]:
                    final_rule = check_transition(s_t4, s_t5, rule2_cand)
                    if final_rule:
                        possible_rule2s.append(((s_t3, s_t4, s_t5), final_rule))
    
    # Rule 3: t=4,5,6
    for s_t4 in size_map[4]:
        for s_t5 in size_map[5]:
            rule3_cand = check_transition(s_t4, s_t5)
            if rule3_cand:
                final_rule = check_transition(s_t5, size_map[6][0], rule3_cand)
                if final_rule:
                    possible_rule3s.append(((s_t4, s_t5, size_map[6][0]), final_rule))

    # 6. Find the unique solution satisfying all conditions
    for r1_seq, rule1 in possible_rule1s:
        for r2_seq, rule2 in possible_rule2s:
            for r3_seq, rule3 in possible_rule3s:
                # Condition: Disjoint States
                all_states = list(r1_seq) + list(r2_seq) + list(r3_seq)
                if len(set(all_states)) == 9:
                    # Condition: Distinct Rules
                    rules_set = set()
                    rules_set.add(tuple(sorted(rule1.items())))
                    rules_set.add(tuple(sorted(rule2.items())))
                    rules_set.add(tuple(sorted(rule3.items())))
                    if len(rules_set) == 3:
                        # 7. Format the output
                        r1_str = "".join(map(str, r1_seq))
                        r2_str = "".join(map(str, r2_seq))
                        r3_str = "".join(map(str, r3_seq))
                        print(f"{{{r1_str},{r2_str},{r3_str}}}")
                        return

solve_ca_puzzle()
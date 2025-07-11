import numpy as np
import itertools

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically searching for consistent
    evolution sequences among the given states.
    """
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

    def parse_state(s):
        lines = s.strip().split('\n')
        return np.array([[int(c) for c in line.split()] for line in lines], dtype=int)

    parsed_states = {i: parse_state(s) for i, s in states_str.items()}
    
    # Memoization cache for transition checks
    transition_cache = {}

    def get_rule_from_transition(grid_t, grid_t1):
        # Align grids by padding the smaller one to match the larger one's dimensions.
        # Since all states are odd-dimensioned and centered, this is straightforward.
        h1, w1 = grid_t.shape
        h2, w2 = grid_t1.shape
        target_dim = max(h1, w1, h2, w2)
        
        pad_h1 = target_dim - h1
        pad_w1 = target_dim - w1
        a_grid_t = np.pad(grid_t, ((pad_h1//2, pad_h1 - pad_h1//2), (pad_w1//2, pad_w1 - pad_w1//2)), 'constant')

        pad_h2 = target_dim - h2
        pad_w2 = target_dim - w2
        a_grid_t1 = np.pad(grid_t1, ((pad_h2//2, pad_h2 - pad_h2//2), (pad_w2//2, pad_w2 - pad_w2//2)), 'constant')

        padded_t = np.pad(a_grid_t, 1, 'constant')
        rule = {}
        
        for r in range(target_dim):
            for c in range(target_dim):
                # Sum the 3x3 Moore neighborhood in the padded t-grid
                neighborhood_sum = np.sum(padded_t[r:r+3, c:c+3])
                outcome = a_grid_t1[r, c]
                
                if neighborhood_sum in rule:
                    if rule[neighborhood_sum] != outcome:
                        return None # Contradiction found
                else:
                    rule[neighborhood_sum] = outcome
        return rule

    def get_rule_memoized(label_t, label_t1):
        if (label_t, label_t1) in transition_cache:
            return transition_cache[(label_t, label_t1)]
        
        grid_t = parsed_states[label_t]
        grid_t1 = parsed_states[label_t1]
        rule = get_rule_from_transition(grid_t, grid_t1)
        transition_cache[(label_t, label_t1)] = rule
        return rule

    def merge_rules(rules):
        merged = {}
        for rule in rules:
            if rule is None: return None
            for s, o in rule.items():
                if s in merged and merged[s] != o:
                    return None # Contradiction during merge
                merged[s] = o
        return merged

    all_state_labels = set(parsed_states.keys())
    # Grid at t=1 is a 3x3 block of 1s. We use label 0 for it.
    parsed_states[0] = np.ones((3, 3), dtype=int)
    
    found_sequences = {}

    # --- Find Rule 1 Sequence (t=2, 3, 4) ---
    for s2, s3, s4 in itertools.permutations(all_state_labels, 3):
        rule_a = get_rule_memoized(0, s2)
        rule_b = get_rule_memoized(s2, s3)
        rule_c = get_rule_memoized(s3, s4)
        if merge_rules([rule_a, rule_b, rule_c]) is not None:
            found_sequences[1] = (s2, s3, s4)
            break
            
    # --- Find Rule 2 Sequence (t=3, 4, 5) ---
    remaining_for_r2 = all_state_labels - set(found_sequences[1])
    for s3, s4, s5 in itertools.permutations(remaining_for_r2, 3):
        rule_a = get_rule_memoized(s3, s4)
        rule_b = get_rule_memoized(s4, s5)
        if merge_rules([rule_a, rule_b]) is not None:
            found_sequences[2] = (s3, s4, s5)
            break

    # --- Find Rule 3 Sequence (t=4, 5, 6) ---
    remaining_for_r3 = remaining_for_r2 - set(found_sequences[2])
    for s4, s5, s6 in itertools.permutations(remaining_for_r3, 3):
        rule_a = get_rule_memoized(s4, s5)
        rule_b = get_rule_memoized(s5, s6)
        if merge_rules([rule_a, rule_b]) is not None:
            found_sequences[3] = (s4, s5, s6)
            break

    r1_labels = "".join(map(str, found_sequences[1]))
    r2_labels = "".join(map(str, found_sequences[2]))
    r3_labels = "".join(map(str, found_sequences[3]))
    
    result_str = f"{{{r1_labels},{r2_labels},{r3_labels}}}"
    print(result_str)

solve_ca_puzzle()
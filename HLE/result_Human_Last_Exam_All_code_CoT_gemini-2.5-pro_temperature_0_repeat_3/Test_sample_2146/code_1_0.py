import numpy as np
from itertools import permutations

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by finding the correct sequence of states for each rule.
    """
    states = {
        1: np.array([[0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]]),
        2: np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]]),
        3: np.array([[0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]]),
        4: np.array([[0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]]),
        5: np.array([[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]]),
        6: np.array([[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]]),
        7: np.array([[1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]]),
        8: np.array([[0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]]),
        9: np.array([[1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]])
    }

    states_by_time = {}
    for label, grid in states.items():
        t = (grid.shape[0] - 1) // 2
        if t not in states_by_time:
            states_by_time[t] = []
        states_by_time[t].append(label)

    def get_rule_from_transition(grid_in, grid_out):
        rule_map = {}
        padded_grid_in = np.pad(grid_in, 2, 'constant')
        
        if grid_out.shape[0] != grid_in.shape[0] + 2:
            return None # Dimension mismatch for standard evolution

        for r in range(grid_out.shape[0]):
            for c in range(grid_out.shape[1]):
                r_in_padded, c_in_padded = r + 1, c + 1
                neighborhood = padded_grid_in[r_in_padded-1:r_in_padded+2, c_in_padded-1:c_in_padded+2]
                s = np.sum(neighborhood)
                next_state = grid_out[r, c]
                if s in rule_map and rule_map[s] != next_state:
                    return None  # Inconsistent rule
                rule_map[s] = next_state
        return rule_map

    def check_sequence_consistency(seq_labels):
        rule_map = {}
        for i in range(len(seq_labels) - 1):
            g_in = states[seq_labels[i]]
            g_out = states[seq_labels[i+1]]
            transition_rule = get_rule_from_transition(g_in, g_out)
            if transition_rule is None:
                return None
            for s, v in transition_rule.items():
                if s in rule_map and rule_map[s] != v:
                    return None
                rule_map[s] = v
        return rule_map

    s3_options = states_by_time[3] # {7, 9}
    s4_options = states_by_time[4] # {4, 6, 8}
    s5_options = states_by_time[5] # {1, 5}

    for p3 in permutations(s3_options):
        s3_r1, s3_r2 = p3
        for p4 in permutations(s4_options):
            s4_r1, s4_r2, s4_r3 = p4
            for p5 in permutations(s5_options):
                s5_r2, s5_r3 = p5
                
                seq1 = (states_by_time[2][0], s3_r1, s4_r1)
                seq2 = (s3_r2, s4_r2, s5_r2)
                seq3 = (s4_r3, s5_r3, states_by_time[6][0])

                rule1 = check_sequence_consistency(seq1)
                if rule1 is None: continue
                
                rule2 = check_sequence_consistency(seq2)
                if rule2 is None: continue

                rule3 = check_sequence_consistency(seq3)
                if rule3 is None: continue

                # Check if rules are distinct
                if rule1 == rule2 or rule1 == rule3 or rule2 == rule3:
                    continue
                
                # Check for conflicts between rules
                all_rules = [rule1, rule2, rule3]
                conflict = False
                for i in range(len(all_rules)):
                    for j in range(i + 1, len(all_rules)):
                        r_a, r_b = all_rules[i], all_rules[j]
                        for s, v in r_a.items():
                            if s in r_b and r_b[s] != v:
                                conflict = True
                                break
                        if conflict: break
                    if conflict: break
                
                if not conflict:
                    r1_labels = "".join(map(str, seq1))
                    r2_labels = "".join(map(str, seq2))
                    r3_labels = "".join(map(str, seq3))
                    print(f"{{{r1_labels},{r2_labels},{r3_labels}}}")
                    return

solve_ca_puzzle()
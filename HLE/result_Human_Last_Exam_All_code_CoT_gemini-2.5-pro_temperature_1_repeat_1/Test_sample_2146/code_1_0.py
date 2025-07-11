import numpy as np
import itertools

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically checking all
    possible assignments of states to rules based on the given constraints.
    """
    states = {
        1: np.array([
            [0,0,0,1,1,1,1,1,0,0,0],
            [0,0,1,0,0,1,0,0,1,0,0],
            [0,1,1,1,1,0,1,1,1,1,0],
            [1,0,1,0,0,1,0,0,1,0,1],
            [1,0,1,0,0,1,0,0,1,0,1],
            [1,1,0,1,1,0,1,1,0,1,1],
            [1,0,1,0,0,1,0,0,1,0,1],
            [1,0,1,0,0,1,0,0,1,0,1],
            [0,1,1,1,1,0,1,1,1,1,0],
            [0,0,1,0,0,1,0,0,1,0,0],
            [0,0,0,1,1,1,1,1,0,0,0]
        ]),
        2: np.array([
            [1,0,0,0,1,0,0,0,1,0,0,0,1],
            [0,1,1,1,0,0,0,0,0,1,1,1,0],
            [0,1,1,1,0,0,0,0,0,1,1,1,0],
            [0,1,1,1,1,1,1,1,1,1,1,1,0],
            [1,0,0,1,1,1,1,1,1,1,0,0,1],
            [0,0,0,1,1,1,1,1,1,1,0,0,0],
            [0,0,0,1,1,1,1,1,1,1,0,0,0],
            [0,0,0,1,1,1,1,1,1,1,0,0,0],
            [1,0,0,1,1,1,1,1,1,1,0,0,1],
            [0,1,1,1,1,1,1,1,1,1,1,1,0],
            [0,1,1,1,0,0,0,0,0,1,1,1,0],
            [0,1,1,1,0,0,0,0,0,1,1,1,0],
            [1,0,0,0,1,0,0,0,1,0,0,0,1]
        ]),
        3: np.array([
            [0,0,1,0,0],
            [0,1,0,1,0],
            [1,0,1,0,1],
            [0,1,0,1,0],
            [0,0,1,0,0]
        ]),
        4: np.array([
            [0,0,0,1,0,1,0,0,0],
            [0,1,1,0,0,0,1,1,0],
            [0,1,1,0,0,0,1,1,0],
            [1,0,0,0,0,0,0,0,1],
            [0,0,0,0,1,0,0,0,0],
            [1,0,0,0,0,0,0,0,1],
            [0,1,1,0,0,0,1,1,0],
            [0,1,1,0,0,0,1,1,0],
            [0,0,0,1,0,1,0,0,0]
        ]),
        5: np.array([
            [1,1,1,0,0,1,0,0,1,1,1],
            [1,0,1,0,1,1,1,0,1,0,1],
            [1,1,1,0,0,1,0,0,1,1,1],
            [0,0,0,1,1,1,1,1,0,0,0],
            [0,1,0,1,1,0,1,1,0,1,0],
            [1,1,1,1,0,0,0,1,1,1,1],
            [0,1,0,1,1,0,1,1,0,1,0],
            [0,0,0,1,1,1,1,1,0,0,0],
            [1,1,1,0,0,1,0,0,1,1,1],
            [1,0,1,0,1,1,1,0,1,0,1],
            [1,1,1,0,0,1,0,0,1,1,1]
        ]),
        6: np.array([
            [1,0,0,0,0,0,0,0,1],
            [0,1,1,0,0,0,1,1,0],
            [0,1,1,1,0,1,1,1,0],
            [0,0,1,1,1,1,1,0,0],
            [0,0,0,1,1,1,0,0,0],
            [0,0,1,1,1,1,1,0,0],
            [0,1,1,1,0,1,1,1,0],
            [0,1,1,0,0,0,1,1,0],
            [1,0,0,0,0,0,0,0,1]
        ]),
        7: np.array([
            [1,1,0,0,0,1,1],
            [1,0,0,0,0,0,1],
            [0,0,1,1,1,0,0],
            [0,0,1,1,1,0,0],
            [0,0,1,1,1,0,0],
            [1,0,0,0,0,0,1],
            [1,1,0,0,0,1,1]
        ]),
        8: np.array([
            [0,0,0,0,1,0,0,0,0],
            [0,0,1,1,0,1,1,0,0],
            [0,1,1,0,0,0,1,1,0],
            [0,1,0,1,1,1,0,1,0],
            [1,0,0,1,0,1,0,0,1],
            [0,1,0,1,1,1,0,1,0],
            [0,1,1,0,0,0,1,1,0],
            [0,0,1,1,0,1,1,0,0],
            [0,0,0,0,1,0,0,0,0]
        ]),
        9: np.array([
            [1,0,1,0,1,0,1],
            [0,1,0,0,0,1,0],
            [1,0,0,1,0,0,1],
            [0,0,1,0,1,0,0],
            [1,0,0,1,0,0,1],
            [0,1,0,0,0,1,0],
            [1,0,1,0,1,0,1]
        ])
    }

    def deduce_rule(g1, g2):
        rule = {}
        # The pattern grows by 1 cell on each side per time step.
        if not (g2.shape[0] == g1.shape[0] + 2 and g2.shape[1] == g1.shape[1] + 2):
            return None 
        
        g1_padded = np.pad(g1, pad_width=1, mode='constant', constant_values=0)
        
        for r in range(g2.shape[0]):
            for c in range(g2.shape[1]):
                neighborhood = g1_padded[r:r+3, c:c+3]
                s = np.sum(neighborhood)
                v = g2[r, c]
                if s in rule and rule[s] != v:
                    return None # Contradiction
                rule[s] = v
        return rule

    def merge_rules(r1, r2):
        merged = r1.copy()
        for s, v in r2.items():
            if s in merged and merged[s] != v:
                return None # Contradiction
            merged[s] = v
        return merged

    # Define the sets of possible labels for each time step
    s3_labels = [7, 9]
    s4_labels = [4, 6, 8]
    s5_labels = [1, 5]

    # Iterate through all permutations of label assignments
    for p_s3 in itertools.permutations(s3_labels):
        for p_s4 in itertools.permutations(s4_labels):
            for p_s5 in itertools.permutations(s5_labels):
                
                # Define the three candidate sequences based on the permutation
                seq1 = (3, p_s3[0], p_s4[0])
                seq2 = (p_s3[1], p_s4[1], p_s5[0])
                seq3 = (p_s4[2], p_s5[1], 2)
                
                # Test Rule 1 sequence
                rule1_part1 = deduce_rule(states[seq1[0]], states[seq1[1]])
                if rule1_part1 is None: continue
                rule1_part2 = deduce_rule(states[seq1[1]], states[seq1[2]])
                if rule1_part2 is None: continue
                full_rule1 = merge_rules(rule1_part1, rule1_part2)
                if full_rule1 is None: continue

                # Test Rule 2 sequence
                rule2_part1 = deduce_rule(states[seq2[0]], states[seq2[1]])
                if rule2_part1 is None: continue
                rule2_part2 = deduce_rule(states[seq2[1]], states[seq2[2]])
                if rule2_part2 is None: continue
                full_rule2 = merge_rules(rule2_part1, rule2_part2)
                if full_rule2 is None: continue
                
                # Test Rule 3 sequence
                rule3_part1 = deduce_rule(states[seq3[0]], states[seq3[1]])
                if rule3_part1 is None: continue
                rule3_part2 = deduce_rule(states[seq3[1]], states[seq3[2]])
                if rule3_part2 is None: continue
                full_rule3 = merge_rules(rule3_part1, rule3_part2)
                if full_rule3 is None: continue

                # Check if all three rules are distinct
                if full_rule1 != full_rule2 and full_rule1 != full_rule3 and full_rule2 != full_rule3:
                    r1_str = "".join(map(str, seq1))
                    r2_str = "".join(map(str, seq2))
                    r3_str = "".join(map(str, seq3))
                    
                    print(f"{{{r1_str},{r2_str},{r3_str}}}")
                    return

solve_ca_puzzle()
import numpy as np
from itertools import permutations

def solve_ca_puzzle():
    """
    This function solves the cellular automaton puzzle by deducing rules from state transitions
    and finding a unique assignment of states to rules that satisfies all conditions.
    """

    # Step 1: Define the nine states as numpy arrays.
    states = {}
    states[1] = np.array([
        [0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],
        [1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],
        [1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],[0,1,1,1,1,0,1,1,1,1,0],
        [0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]
    ])
    states[2] = np.array([
        [1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],
        [0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],
        [0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],
        [0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],
        [1,0,0,0,1,0,0,0,1,0,0,0,1]
    ])
    states[3] = np.array([
        [0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]
    ])
    states[4] = np.array([
        [0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],
        [1,0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],
        [0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[0,0,0,1,0,1,0,0,0]
    ])
    states[5] = np.array([
        [1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],
        [0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],
        [0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],
        [1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]
    ])
    states[6] = np.array([
        [1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],
        [0,0,1,1,1,1,1,0,0],[0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],
        [0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1]
    ])
    states[7] = np.array([
        [1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],
        [0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]
    ])
    states[8] = np.array([
        [0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],
        [0,1,0,1,1,1,0,1,0],[1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],
        [0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],[0,0,0,0,1,0,0,0,0]
    ])
    states[9] = np.array([
        [1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],
        [1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]
    ])

    # Step 2: Define helper functions to deduce and merge rules.
    def deduce_rule(grid_t, grid_t_plus_1):
        h_t, w_t = grid_t.shape
        h_tp1, w_tp1 = grid_t_plus_1.shape
        # A pattern at time t is size (2t+1). At t+1, it's (2(t+1)+1) = 2t+3.
        # The dimension difference should be 2.
        if h_tp1 != h_t + 2 or w_tp1 != w_t + 2: return None
        
        grid_t_padded = np.pad(grid_t, ((1, 1), (1, 1)), 'constant')
        rule = {}
        for r in range(h_tp1):
            for c in range(w_tp1):
                s = np.sum(grid_t_padded[r:r+3, c:c+3])
                next_state = grid_t_plus_1[r, c]
                if s in rule and rule[s] != next_state: return None
                rule[s] = next_state
        return rule

    def merge_rules(rule1, rule2):
        merged = rule1.copy()
        for s, state in rule2.items():
            if s in merged and merged[s] != state: return None
            merged[s] = state
        return merged

    # Step 3: Find all valid one-step transitions and two-step chains.
    transitions = {}
    for i in range(1, 10):
        for j in range(1, 10):
            # Check if dimensions correspond to a one-step evolution
            if states[i].shape[0] + 2 == states[j].shape[0]:
                rule = deduce_rule(states[i], states[j])
                if rule is not None:
                    transitions[(i, j)] = rule
    
    chains = {}
    for (s1, s2), rule1 in transitions.items():
        for (s2_prime, s3), rule2 in transitions.items():
            if s2 == s2_prime:
                merged = merge_rules(rule1, rule2)
                if merged is not None:
                    chains[(s1, s2, s3)] = merged

    # Step 4: Iterate through permutations of state assignments to find the unique solution.
    # Based on state sizes: t=2 -> 5x5 (#3), t=3 -> 7x7 (#7,#9), t=4 -> 9x9 (#4,#6,#8),
    # t=5 -> 11x11 (#1,#5), t=6 -> 13x13 (#2).
    T3_states, T4_states, T5_states = [7, 9], [4, 6, 8], [1, 5]
    
    solution = None
    for s3_R1, s3_R2 in permutations(T3_states):
        for s4_R1, s4_R2, s4_R3 in permutations(T4_states):
            for s5_R2, s5_R3 in permutations(T5_states):
                # Rule 1: t=2,3,4. Rule 2: t=3,4,5. Rule 3: t=4,5,6
                chain1 = (3, s3_R1, s4_R1)
                chain2 = (s3_R2, s4_R2, s5_R2)
                chain3 = (s4_R3, s5_R3, 2)
                
                if chain1 in chains and chain2 in chains and chain3 in chains:
                    rule1, rule2, rule3 = chains[chain1], chains[chain2], chains[chain3]
                    # Check for distinct rules
                    if rule1 != rule2 and rule1 != rule3 and rule2 != rule3:
                        solution = (chain1, chain2, chain3)
                        break
            if solution: break
        if solution: break

    # Step 5: Format and print the final answer.
    R1_Labels = "".join(map(str, solution[0]))
    R2_Labels = "".join(map(str, solution[1]))
    R3_Labels = "".join(map(str, solution[2]))
    
    print(f"{{{R1_Labels},{R2_Labels},{R3_Labels}}}")

solve_ca_puzzle()
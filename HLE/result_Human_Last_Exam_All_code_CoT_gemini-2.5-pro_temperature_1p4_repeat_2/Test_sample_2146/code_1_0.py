import numpy as np
import itertools

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by identifying the state sequences for three distinct rules.
    """

    # 1. Represent the states as numpy arrays.
    states = {
        1: np.array([
            [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0], [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
            [0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0], [1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1],
            [1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1], [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1],
            [1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1], [1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1],
            [0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0], [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
            [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
        ]),
        2: np.array([
            [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0],
            [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0], [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
            [1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1], [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
            [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1], [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
            [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0], [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0],
            [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
        ]),
        3: np.array([
            [0, 0, 1, 0, 0], [0, 1, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [0, 0, 1, 0, 0],
        ]),
        4: np.array([
            [0, 0, 0, 1, 0, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0, 1, 1, 0],
            [0, 1, 1, 0, 0, 0, 1, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 1, 1, 0, 0, 0, 1, 1, 0], [0, 1, 1, 0, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 1, 0, 0, 0],
        ]),
        5: np.array([
            [1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1], [1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1],
            [1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
            [0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0], [1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1],
            [0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0], [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
            [1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1], [1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1],
            [1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1],
        ]),
        6: np.array([
            [1, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 1, 0, 0, 0, 1, 1, 0],
            [0, 1, 1, 1, 0, 1, 1, 1, 0], [0, 0, 1, 1, 1, 1, 1, 0, 0],
            [0, 0, 0, 1, 1, 1, 0, 0, 0], [0, 0, 1, 1, 1, 1, 1, 0, 0],
            [0, 1, 1, 1, 0, 1, 1, 1, 0], [0, 1, 1, 0, 0, 0, 1, 1, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 1],
        ]),
        7: np.array([
            [1, 1, 0, 0, 0, 1, 1], [1, 0, 0, 0, 0, 0, 1], [0, 0, 1, 1, 1, 0, 0],
            [0, 0, 1, 1, 1, 0, 0], [0, 0, 1, 1, 1, 0, 0], [1, 0, 0, 0, 0, 0, 1],
            [1, 1, 0, 0, 0, 1, 1],
        ]),
        8: np.array([
            [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 1, 0, 1, 1, 0, 0],
            [0, 1, 1, 0, 0, 0, 1, 1, 0], [0, 1, 0, 1, 1, 1, 0, 1, 0],
            [1, 0, 0, 1, 0, 1, 0, 0, 1], [0, 1, 0, 1, 1, 1, 0, 1, 0],
            [0, 1, 1, 0, 0, 0, 1, 1, 0], [0, 0, 1, 1, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
        ]),
        9: np.array([
            [1, 0, 1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1, 0], [1, 0, 0, 1, 0, 0, 1],
            [0, 0, 1, 0, 1, 0, 0], [1, 0, 0, 1, 0, 0, 1], [0, 1, 0, 0, 0, 1, 0],
            [1, 0, 1, 0, 1, 0, 1],
        ]),
    }

    # 2 & 3. Group states by time step based on grid size (2t+1).
    states_by_time = {
        2: [s for s, g in states.items() if g.shape[0] == 5],  # 2*2+1
        3: [s for s, g in states.items() if g.shape[0] == 7],  # 2*3+1
        4: [s for s, g in states.items() if g.shape[0] == 9],  # 2*4+1
        5: [s for s, g in states.items() if g.shape[0] == 11], # 2*5+1
        6: [s for s, g in states.items() if g.shape[0] == 13], # 2*6+1
    }

    # 4. Helper functions for rule deduction and verification.
    def get_rule_from_pair(grid_before, grid_after):
        if grid_after.shape[0] != grid_before.shape[0] + 2: return None
        padded_before = np.pad(grid_before, pad_width=2, mode='constant', constant_values=0)
        rule = {}
        for r, c in np.ndindex(grid_after.shape):
            neighborhood_sum = int(np.sum(padded_before[r:r+3, c:c+3]))
            next_state = grid_after[r, c]
            if neighborhood_sum in rule and rule[neighborhood_sum] != next_state: return None
            rule[neighborhood_sum] = next_state
        return rule

    def verify_evolution(grid_before, grid_after, rule):
        if grid_after.shape[0] != grid_before.shape[0] + 2: return False
        padded_before = np.pad(grid_before, pad_width=2, mode='constant', constant_values=0)
        for r, c in np.ndindex(grid_after.shape):
            neighborhood_sum = int(np.sum(padded_before[r:r+3, c:c+3]))
            if grid_after[r, c] != rule.get(neighborhood_sum, 0): return False
        return True

    # 5. Systematic search for rule sequences.
    found_sequences = {}
    all_states = set(states.keys())

    # Find Rule 3 (t=4, 5, 6) first as it's most constrained.
    for s4, s5 in itertools.product(states_by_time[4], states_by_time[5]):
        rule = get_rule_from_pair(states[s4], states[s5])
        if rule and verify_evolution(states[s5], states[states_by_time[6][0]], rule):
            found_sequences[3] = (s4, s5, states_by_time[6][0])
            break
    
    # Find Rule 1 (t=2, 3, 4) from remaining states.
    used = set(found_sequences[3])
    rem_s3 = [s for s in states_by_time[3] if s not in used]
    rem_s4 = [s for s in states_by_time[4] if s not in used]
    s2 = states_by_time[2][0]
    for s3, s4 in itertools.product(rem_s3, rem_s4):
        rule = get_rule_from_pair(states[s2], states[s3])
        if rule and verify_evolution(states[s3], states[s4], rule):
            found_sequences[1] = (s2, s3, s4)
            break

    # Rule 2 must be the last remaining states.
    used.update(found_sequences[1])
    remaining = tuple(sorted(list(all_states - used)))
    found_sequences[2] = remaining
    
    # 6. Assemble and print the final answer string.
    r1_labels = "".join(map(str, found_sequences[1]))
    r2_labels = "".join(map(str, found_sequences[2]))
    r3_labels = "".join(map(str, found_sequences[3]))

    print(f"{{{r1_labels},{r2_labels},{r3_labels}}}")

solve_ca_puzzle()
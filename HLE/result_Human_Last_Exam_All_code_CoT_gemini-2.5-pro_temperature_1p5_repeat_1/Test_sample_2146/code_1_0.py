import numpy as np
from collections import defaultdict
import itertools

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically identifying the rules and state sequences.
    """
    
    # Step 1: Define the states from the problem description
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
    
    # Step 2: Categorize states by time step based on grid size
    states_by_time = defaultdict(list)
    for i, s in states.items():
        time = (s.shape[0] - 1) // 2
        states_by_time[time].append(i)

    def get_rule_from_transitions(sequence):
        """Checks if a sequence of states is valid and returns the implied rule."""
        on_sums = set()
        off_sums = set()
        
        for i in range(len(sequence) - 1):
            prev_state = states[sequence[i]]
            next_state = states[sequence[i+1]]

            # Dimension check
            if next_state.shape[0] != prev_state.shape[0] + 2:
                return None
            
            padded_prev = np.pad(prev_state, pad_width=1, mode='constant', constant_values=0)
            
            for r in range(next_state.shape[0]):
                for c in range(next_state.shape[1]):
                    s = np.sum(padded_prev[r:r+3, c:c+3])
                    next_val = next_state[r, c]
                    
                    if next_val == 1:
                        if s in off_sums: return None  # Contradiction
                        on_sums.add(s)
                    else:
                        if s in on_sums: return None  # Contradiction
                        off_sums.add(s)
        return (on_sums, off_sums)

    solutions = {}
    
    # Step 3 & 4: Search for rule sequences
    
    # Find Rule 3 (t=4,5,6)
    r3_candidates = itertools.product(states_by_time[4], states_by_time[5], states_by_time[6])
    for cand in r3_candidates:
        if get_rule_from_transitions(cand) is not None:
            solutions[3] = cand
            break
            
    # Find Rule 1 (t=2,3,4) from remaining states
    used_states = set(solutions[3])
    r1_t3_pool = [s for s in states_by_time[3] if s not in used_states]
    r1_t4_pool = [s for s in states_by_time[4] if s not in used_states]
    r1_candidates = itertools.product(states_by_time[2], r1_t3_pool, r1_t4_pool)
    for cand in r1_candidates:
        if get_rule_from_transitions(cand) is not None:
            solutions[1] = cand
            break

    # Find Rule 2 (t=3,4,5) from the rest
    used_states.update(solutions[1])
    r2_t3_pool = [s for s in states_by_time[3] if s not in used_states]
    r2_t4_pool = [s for s in states_by_time[4] if s not in used_states]
    r2_t5_pool = [s for s in states_by_time[5] if s not in used_states]
    solutions[2] = (r2_t3_pool[0], r2_t4_pool[0], r2_t5_pool[0])

    # Step 5: Format the final answer
    r1_labels = "".join(map(str, solutions[1]))
    r2_labels = "".join(map(str, solutions[2]))
    r3_labels = "".join(map(str, solutions[3]))
    
    final_answer_string = f"{{{r1_labels},{r2_labels},{r3_labels}}}"
    print(final_answer_string)

solve_ca_puzzle()
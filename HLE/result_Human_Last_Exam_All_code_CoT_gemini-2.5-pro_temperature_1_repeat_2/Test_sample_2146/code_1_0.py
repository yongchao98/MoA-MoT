import numpy as np

def solve_cellular_automaton():
    """
    This function solves the cellular automaton puzzle by:
    1.  Representing the 9 given states as numpy arrays.
    2.  Brute-forcing all 512 possible 9-neighbor totalistic rules where an all-zero neighborhood remains zero.
    3.  For each rule, simulating the evolution from a single '1' cell up to t=6.
    4.  Collecting all rules that produce a sequence of states matching the time-step constraints for Rule 1, 2, or 3.
    5.  Searching through the collected matches to find a unique combination of three distinct rules that uses each of the 9 state labels exactly once.
    6.  Formatting and printing the final answer.
    """
    # Step 1: Define the 9 states as numpy arrays
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

    # Group states by time step (determined by dimension: t = (dim-1)/2)
    # The dictionary keys are matrix data (as bytes), values are state labels
    states_by_time_bytes = {
        2: {m.tobytes(): l for l, m in {3: states[3]}.items()},
        3: {m.tobytes(): l for l, m in {7: states[7], 9: states[9]}.items()},
        4: {m.tobytes(): l for l, m in {4: states[4], 6: states[6], 8: states[8]}.items()},
        5: {m.tobytes(): l for l, m in {1: states[1], 5: states[5]}.items()},
        6: {m.tobytes(): l for l, m in {2: states[2]}.items()}
    }

    # CA simulation utility functions
    def crop(grid):
        if not np.any(grid): return np.array([[0]])
        rows, cols = np.where(grid >= 1)
        if rows.size == 0: return np.array([[0]])
        return grid[rows.min():rows.max()+1, cols.min():cols.max()+1]

    def simulate_step(grid, rule):
        padded_grid = np.pad(grid, pad_width=2, mode='constant', constant_values=0)
        new_grid = np.zeros_like(padded_grid)
        for r in range(1, padded_grid.shape[0] - 1):
            for c in range(1, padded_grid.shape[1] - 1):
                neighborhood_sum = np.sum(padded_grid[r-1:r+2, c-1:c+2])
                if neighborhood_sum < 10:
                    new_grid[r, c] = rule[int(neighborhood_sum)]
        return crop(new_grid)

    # Find all rules that produce valid sequences
    all_matches = {1: [], 2: [], 3: []}
    for i in range(512):  # 2^9 possibilities for sums 1-9
        rule_str = format(i, '09b')
        rule = tuple([0] + [int(bit) for bit in reversed(rule_str)]) # rule[0]=0
        
        grid = np.array([[1]])
        gen_states_bytes = {}
        for t in range(1, 7):
            grid = simulate_step(grid, rule)
            gen_states_bytes[t] = grid.tobytes()

        # Check for Rule 1 pattern: t=2,3,4
        if (gen_states_bytes[2] in states_by_time_bytes[2] and
            gen_states_bytes[3] in states_by_time_bytes[3] and
            gen_states_bytes[4] in states_by_time_bytes[4]):
            labels = (states_by_time_bytes[2][gen_states_bytes[2]],
                      states_by_time_bytes[3][gen_states_bytes[3]],
                      states_by_time_bytes[4][gen_states_bytes[4]])
            all_matches[1].append({'rule': rule, 'labels': labels})

        # Check for Rule 2 pattern: t=3,4,5
        if (gen_states_bytes[3] in states_by_time_bytes[3] and
            gen_states_bytes[4] in states_by_time_bytes[4] and
            gen_states_bytes[5] in states_by_time_bytes[5]):
            labels = (states_by_time_bytes[3][gen_states_bytes[3]],
                      states_by_time_bytes[4][gen_states_bytes[4]],
                      states_by_time_bytes[5][gen_states_bytes[5]])
            all_matches[2].append({'rule': rule, 'labels': labels})

        # Check for Rule 3 pattern: t=4,5,6
        if (gen_states_bytes[4] in states_by_time_bytes[4] and
            gen_states_bytes[5] in states_by_time_bytes[5] and
            gen_states_bytes[6] in states_by_time_bytes[6]):
            labels = (states_by_time_bytes[4][gen_states_bytes[4]],
                      states_by_time_bytes[5][gen_states_bytes[5]],
                      states_by_time_bytes[6][gen_states_bytes[6]])
            all_matches[3].append({'rule': rule, 'labels': labels})
    
    # Find the unique combination of 3 distinct rules that uses all 9 labels
    solution = {}
    for m1 in all_matches[1]:
        for m2 in all_matches[2]:
            for m3 in all_matches[3]:
                # Condition 1: Rules must be distinct
                if len({m1['rule'], m2['rule'], m3['rule']}) == 3:
                    # Condition 2: All 9 labels must be used exactly once
                    all_labels = list(m1['labels']) + list(m2['labels']) + list(m3['labels'])
                    if len(set(all_labels)) == 9:
                        solution[1] = "".join(map(str, m1['labels']))
                        solution[2] = "".join(map(str, m2['labels']))
                        solution[3] = "".join(map(str, m3['labels']))
                        break
            if solution: break
        if solution: break
        
    # Format and print the final result
    result_string = f"{{{solution[1]},{solution[2]},{solution[3]}}}"
    print(result_string)

solve_cellular_automaton()
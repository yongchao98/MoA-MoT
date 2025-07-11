import numpy as np

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically finding the state sequences
    for the three distinct rules.
    """
    states_raw = {
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

    states = {}
    for i, s in states_raw.items():
        lines = s.strip().split('\n')
        grid = [list(map(int, line.split())) for line in lines]
        states[i] = np.array(grid, dtype=int)

    states_by_time = {}
    for id, grid in states.items():
        t = (grid.shape[0] - 1) // 2
        states_by_time.setdefault(t, []).append(id)

    memo_sums = {}
    def get_sums(prev_grid, next_grid):
        key = (prev_grid.tobytes(), next_grid.tobytes())
        if key in memo_sums:
            return memo_sums[key]

        padded_prev = np.pad(prev_grid, 2, 'constant')
        sums_for_1 = set()
        sums_for_0 = set()
        h, w = next_grid.shape
        for r in range(h):
            for c in range(w):
                neighborhood_sum = np.sum(padded_prev[r:r+3, c:c+3])
                if next_grid[r, c] == 1:
                    sums_for_1.add(neighborhood_sum)
                else:
                    sums_for_0.add(neighborhood_sum)
        
        memo_sums[key] = (sums_for_1, sums_for_0)
        return sums_for_1, sums_for_0

    # --- Find Rule 1 (t=2, 3, 4) ---
    rule1_seq, rule1_constraints = None, None
    g2_id = states_by_time[2][0]
    for g3_id in states_by_time[3]:
        s1a, s0a = get_sums(states[g2_id], states[g3_id])
        if not s1a.isdisjoint(s0a): continue
        for g4_id in states_by_time[4]:
            s1b, s0b = get_sums(states[g3_id], states[g4_id])
            s1_total, s0_total = s1a.union(s1b), s0a.union(s0b)
            if s1_total.isdisjoint(s0_total):
                rule1_seq = (g2_id, g3_id, g4_id)
                rule1_constraints = (s1_total, s0_total)
                break
        if rule1_seq: break
    
    # --- Find Rule 2 (t=3, 4, 5) ---
    used_ids = set(rule1_seq)
    rem_t3 = [i for i in states_by_time[3] if i not in used_ids]
    rem_t4 = [i for i in states_by_time[4] if i not in used_ids]
    rem_t5 = states_by_time[5]
    
    rule2_seq, rule2_constraints = None, None
    g3_id = rem_t3[0]
    for g4_id in rem_t4:
        s1a, s0a = get_sums(states[g3_id], states[g4_id])
        if not s1a.isdisjoint(s0a): continue
        for g5_id in rem_t5:
            s1b, s0b = get_sums(states[g4_id], states[g5_id])
            s1_total, s0_total = s1a.union(s1b), s0a.union(s0b)
            if s1_total.isdisjoint(s0_total) and (s1_total, s0_total) != rule1_constraints:
                rule2_seq = (g3_id, g4_id, g5_id)
                rule2_constraints = (s1_total, s0_total)
                break
        if rule2_seq: break

    # --- Find Rule 3 (t=4, 5, 6) ---
    used_ids.update(rule2_seq)
    g4_id = [i for i in states_by_time[4] if i not in used_ids][0]
    g5_id = [i for i in states_by_time[5] if i not in used_ids][0]
    g6_id = states_by_time[6][0]
    rule3_seq = (g4_id, g5_id, g6_id)
    
    # Format and print the final answer
    r1_str = "".join(map(str, rule1_seq))
    r2_str = "".join(map(str, rule2_seq))
    r3_str = "".join(map(str, rule3_seq))

    print(f"Rule 1 corresponds to states at t=2,3,4 with labels: {' -> '.join(map(str,rule1_seq))}")
    print(f"Rule 2 corresponds to states at t=3,4,5 with labels: {' -> '.join(map(str,rule2_seq))}")
    print(f"Rule 3 corresponds to states at t=4,5,6 with labels: {' -> '.join(map(str,rule3_seq))}")
    print(f"\nFinal answer string:")
    print(f"{{{r1_str},{r2_str},{r3_str}}}")

solve_ca_puzzle()
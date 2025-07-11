import numpy as np
import itertools

def solve():
    """
    This function solves the cellular automaton puzzle by identifying the rule 
    and state sequence for three different rules. The logic follows the
    step-by-step plan outlined above.
    """
    states_data = {
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

    def derive_rule(prev_state, next_state):
        live_sums, dead_sums = set(), set()
        nh, nw = next_state.shape
        ph, pw = prev_state.shape
        padded_prev = np.zeros((nh + 2, nw + 2), dtype=int)
        start_h, start_w = (nh + 2 - ph) // 2, (nw + 2 - pw) // 2
        padded_prev[start_h:start_h + ph, start_w:start_w + pw] = prev_state
        for r, c in np.ndindex(next_state.shape):
            s = np.sum(padded_prev[r:r + 3, c:c + 3])
            if next_state[r, c] == 1: live_sums.add(s)
            else: dead_sums.add(s)
        return (live_sums, dead_sums) if live_sums.isdisjoint(dead_sums) else None

    def combine_rules(rule1, rule2):
        live1, dead1 = rule1
        live2, dead2 = rule2
        if not live1.isdisjoint(dead2) or not live2.isdisjoint(dead1): return None
        return (live1 | live2, dead1 | dead2)

    t_map = {2: [3], 3: [7, 9], 4: [4, 6, 8], 5: [1, 5], 6: [2]}
    rule_cache = {}

    def get_rule(s1_id, s2_id):
        if (s1_id, s2_id) not in rule_cache:
            rule_cache[(s1_id, s2_id)] = derive_rule(states_data[s1_id], states_data[s2_id])
        return rule_cache[(s1_id, s2_id)]

    for t3_s1, t3_s2 in itertools.permutations(t_map[3]):
        for t4_s1, t4_s2, t4_s3 in itertools.permutations(t_map[4]):
            for t5_s2, t5_s3 in itertools.permutations(t_map[5]):
                path1, path2, path3 = (t_map[2][0], t3_s1, t4_s1), (t3_s2, t4_s2, t5_s2), (t4_s3, t5_s3, t_map[6][0])
                
                r1p1, r1p2 = get_rule(path1[0], path1[1]), get_rule(path1[1], path1[2])
                if not r1p1 or not r1p2: continue
                full_rule1 = combine_rules(r1p1, r1p2)
                if not full_rule1: continue

                r2p1, r2p2 = get_rule(path2[0], path2[1]), get_rule(path2[1], path2[2])
                if not r2p1 or not r2p2: continue
                full_rule2 = combine_rules(r2p1, r2p2)
                if not full_rule2: continue

                r3p1, r3p2 = get_rule(path3[0], path3[1]), get_rule(path3[1], path3[2])
                if not r3p1 or not r3p2: continue
                full_rule3 = combine_rules(r3p1, r3p2)
                if not full_rule3: continue

                if len({frozenset(r[0]) for r in [full_rule1, full_rule2, full_rule3]}) == 3:
                    r1_labels = "".join(map(str, path1))
                    r2_labels = "".join(map(str, path2))
                    r3_labels = "".join(map(str, path3))
                    print(f"{{{r1_labels},{r2_labels},{r3_labels}}}")
                    return

solve()
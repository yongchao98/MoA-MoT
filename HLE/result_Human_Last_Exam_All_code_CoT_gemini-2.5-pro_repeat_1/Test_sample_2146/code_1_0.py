import numpy as np
from itertools import permutations

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by identifying the sequence of states for three different rules.
    """
    states = {
        1: np.array([
            [0,0,0,1,1,1,1,1,0,0,0],[0,0,1,0,0,1,0,0,1,0,0],[0,1,1,1,1,0,1,1,1,1,0],[1,0,1,0,0,1,0,0,1,0,1],
            [1,0,1,0,0,1,0,0,1,0,1],[1,1,0,1,1,0,1,1,0,1,1],[1,0,1,0,0,1,0,0,1,0,1],[1,0,1,0,0,1,0,0,1,0,1],
            [0,1,1,1,1,0,1,1,1,1,0],[0,0,1,0,0,1,0,0,1,0,0],[0,0,0,1,1,1,1,1,0,0,0]
        ], dtype=int),
        2: np.array([
            [1,0,0,0,1,0,0,0,1,0,0,0,1],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],
            [0,1,1,1,1,1,1,1,1,1,1,1,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],[0,0,0,1,1,1,1,1,1,1,0,0,0],
            [0,0,0,1,1,1,1,1,1,1,0,0,0],[0,0,0,1,1,1,1,1,1,1,0,0,0],[1,0,0,1,1,1,1,1,1,1,0,0,1],
            [0,1,1,1,1,1,1,1,1,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0,0,1,1,1,0],
            [1,0,0,0,1,0,0,0,1,0,0,0,1]
        ], dtype=int),
        3: np.array([
            [0,0,1,0,0],[0,1,0,1,0],[1,0,1,0,1],[0,1,0,1,0],[0,0,1,0,0]
        ], dtype=int),
        4: np.array([
            [0,0,0,1,0,1,0,0,0],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],[1,0,0,0,0,0,0,0,1],
            [0,0,0,0,1,0,0,0,0],[1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,0,0,0,1,1,0],
            [0,0,0,1,0,1,0,0,0]
        ], dtype=int),
        5: np.array([
            [1,1,1,0,0,1,0,0,1,1,1],[1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1],
            [0,0,0,1,1,1,1,1,0,0,0],[0,1,0,1,1,0,1,1,0,1,0],[1,1,1,1,0,0,0,1,1,1,1],
            [0,1,0,1,1,0,1,1,0,1,0],[0,0,0,1,1,1,1,1,0,0,0],[1,1,1,0,0,1,0,0,1,1,1],
            [1,0,1,0,1,1,1,0,1,0,1],[1,1,1,0,0,1,0,0,1,1,1]
        ], dtype=int),
        6: np.array([
            [1,0,0,0,0,0,0,0,1],[0,1,1,0,0,0,1,1,0],[0,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,1,0,0],
            [0,0,0,1,1,1,0,0,0],[0,0,1,1,1,1,1,0,0],[0,1,1,1,0,1,1,1,0],[0,1,1,0,0,0,1,1,0],
            [1,0,0,0,0,0,0,0,1]
        ], dtype=int),
        7: np.array([
            [1,1,0,0,0,1,1],[1,0,0,0,0,0,1],[0,0,1,1,1,0,0],[0,0,1,1,1,0,0],
            [0,0,1,1,1,0,0],[1,0,0,0,0,0,1],[1,1,0,0,0,1,1]
        ], dtype=int),
        8: np.array([
            [0,0,0,0,1,0,0,0,0],[0,0,1,1,0,1,1,0,0],[0,1,1,0,0,0,1,1,0],[0,1,0,1,1,1,0,1,0],
            [1,0,0,1,0,1,0,0,1],[0,1,0,1,1,1,0,1,0],[0,1,1,0,0,0,1,1,0],[0,0,1,1,0,1,1,0,0],
            [0,0,0,0,1,0,0,0,0]
        ], dtype=int),
        9: np.array([
            [1,0,1,0,1,0,1],[0,1,0,0,0,1,0],[1,0,0,1,0,0,1],[0,0,1,0,1,0,0],
            [1,0,0,1,0,0,1],[0,1,0,0,0,1,0],[1,0,1,0,1,0,1]
        ], dtype=int),
    }

    def find_rule(state_before, state_after):
        live_sums, dead_sums = set(), set()
        padded_before = np.pad(state_before, 1, 'constant')
        if padded_before.shape != state_after.shape: return None

        rows, cols = state_after.shape
        for r in range(rows):
            for c in range(cols):
                s = np.sum(padded_before[r:r+3, c:c+3])
                if state_after[r, c] == 1:
                    if s in dead_sums: return None
                    live_sums.add(s)
                else:
                    if s in live_sums: return None
                    dead_sums.add(s)
        return {'live': live_sums, 'dead': dead_sums}

    def combine_rules(rule1, rule2):
        if not rule1 or not rule2: return None
        live_sums = rule1['live'].union(rule2['live'])
        dead_sums = rule1['dead'].union(rule2['dead'])
        if not live_sums.isdisjoint(dead_sums): return None
        return {'live': live_sums, 'dead': dead_sums}

    state_groups = {t: [] for t in range(2, 7)}
    for label, matrix in states.items():
        t = (matrix.shape[0] - 1) // 2
        if t in state_groups:
            state_groups[t].append(label)

    rule_cache = {}
    for t in range(2, 6):
        for s_before in state_groups[t]:
            for s_after in state_groups[t+1]:
                rule = find_rule(states[s_before], states[s_after])
                if rule: rule_cache[(s_before, s_after)] = rule

    def solve():
        s2_pool, s3_pool, s4_pool, s5_pool, s6_pool = (state_groups[2], state_groups[3], state_groups[4], state_groups[5], state_groups[6])
        
        # Rule 1: t=2,3,4. Starts with state #3
        s2_r1 = 3
        for s3_r1 in s3_pool:
            for s4_r1 in s4_pool:
                r1_seq = (s2_r1, s3_r1, s4_r1)
                rule1_trans1 = rule_cache.get((s2_r1, s3_r1))
                rule1_trans2 = rule_cache.get((s3_r1, s4_r1))
                if not combine_rules(rule1_trans1, rule1_trans2): continue
                
                # Rule 2: t=3,4,5
                used1 = set(r1_seq)
                for s3_r2 in [s for s in s3_pool if s not in used1]:
                    for s4_r2 in [s for s in s4_pool if s not in used1]:
                        for s5_r2 in [s for s in s5_pool if s not in used1]:
                            r2_seq = (s3_r2, s4_r2, s5_r2)
                            rule2_trans1 = rule_cache.get((s3_r2, s4_r2))
                            rule2_trans2 = rule_cache.get((s4_r2, s5_r2))
                            if not combine_rules(rule2_trans1, rule2_trans2): continue
                            
                            # Rule 3: t=4,5,6
                            used2 = used1.union(r2_seq)
                            s4_r3 = next((s for s in s4_pool if s not in used2), None)
                            s5_r3 = next((s for s in s5_pool if s not in used2), None)
                            s6_r3 = next((s for s in s6_pool if s not in used2), None)
                            
                            if None in [s4_r3, s5_r3, s6_r3]: continue

                            r3_seq = (s4_r3, s5_r3, s6_r3)
                            rule3_trans1 = rule_cache.get((s4_r3, s5_r3))
                            rule3_trans2 = rule_cache.get((s5_r3, s6_r3))
                            if combine_rules(rule3_trans1, rule3_trans2):
                                return r1_seq, r2_seq, r3_seq
        return None, None, None

    r1_sol, r2_sol, r3_sol = solve()
    
    r1_str = "".join(map(str, r1_sol))
    r2_str = "".join(map(str, r2_sol))
    r3_str = "".join(map(str, r3_sol))
    
    print(f"{{{r1_str},{r2_str},{r3_str}}}")

solve_ca_puzzle()
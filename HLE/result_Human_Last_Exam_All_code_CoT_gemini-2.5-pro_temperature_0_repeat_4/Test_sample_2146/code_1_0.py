import itertools

def solve_ca_puzzle():
    """
    Solves the cellular automaton puzzle by systematically searching for consistent state sequences.
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
    states = {i: [[int(c) for c in line.split()] for line in s.strip().split('\n')] for i, s in states_raw.items()}

    def get_rule_from_transition(grid1, grid2):
        live_sums = {0} 
        dead_sums = set()
        h1, w1 = len(grid1), len(grid1[0])
        h2, w2 = len(grid2), len(grid2[0])
        if h2 != h1 + 2 or w2 != w1 + 2: return None

        padded_grid1 = [[1] * (w1 + 2) for _ in range(h1 + 2)]
        for r in range(h1):
            for c in range(w1):
                padded_grid1[r + 1][c + 1] = grid1[r][c]

        for r in range(h2):
            for c in range(w2):
                current_sum = 0
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        nr, nc = r + dr, c + dc
                        current_sum += padded_grid1[nr][nc] if 0 <= nr < h2 and 0 <= nc < w2 else 1
                
                if grid2[r][c] == 1:
                    if current_sum in dead_sums: return None
                    live_sums.add(current_sum)
                else:
                    if current_sum in live_sums: return None
                    dead_sums.add(current_sum)
        return (frozenset(live_sums), frozenset(dead_sums))

    def combine_rules(rule1, rule2):
        if not rule1 or not rule2: return None
        live1, dead1 = rule1
        live2, dead2 = rule2
        if live1.intersection(dead2) or live2.intersection(dead1): return None
        return (live1.union(live2), dead1.union(dead2))

    states_by_time = {3: [7, 9], 4: [4, 6, 8], 5: [1, 5]}
    s3_perms = list(itertools.permutations(states_by_time[3]))
    s4_perms = list(itertools.permutations(states_by_time[4]))
    s5_perms = list(itertools.permutations(states_by_time[5]))

    for s3_p in s3_perms:
        for s4_p in s4_perms:
            for s5_p in s5_perms:
                rule1_seq = (3, s3_p[0], s4_p[0])
                rule2_seq = (s3_p[1], s4_p[1], s5_p[0])
                rule3_seq = (s4_p[2], s5_p[1], 2)

                r1_c1 = get_rule_from_transition(states[rule1_seq[0]], states[rule1_seq[1]])
                r1_c2 = get_rule_from_transition(states[rule1_seq[1]], states[rule1_seq[2]])
                rule1_final = combine_rules(r1_c1, r1_c2)
                if rule1_final is None: continue

                r2_c1 = get_rule_from_transition(states[rule2_seq[0]], states[rule2_seq[1]])
                r2_c2 = get_rule_from_transition(states[rule2_seq[1]], states[rule2_seq[2]])
                rule2_final = combine_rules(r2_c1, r2_c2)
                if rule2_final is None: continue

                r3_c1 = get_rule_from_transition(states[rule3_seq[0]], states[rule3_seq[1]])
                r3_c2 = get_rule_from_transition(states[rule3_seq[1]], states[rule3_seq[2]])
                rule3_final = combine_rules(r3_c1, r3_c2)
                if rule3_final is None: continue
                
                if len({rule1_final, rule2_final, rule3_final}) == 3:
                    r1_str = "".join(map(str, rule1_seq))
                    r2_str = "".join(map(str, rule2_seq))
                    r3_str = "".join(map(str, rule3_seq))
                    print(f"{{{r1_str},{r2_str},{r3_str}}}")
                    return

solve_ca_puzzle()